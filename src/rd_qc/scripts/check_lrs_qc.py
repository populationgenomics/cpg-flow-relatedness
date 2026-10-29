"""
Checks long-read CRAM QC metrics in MultiQC output against per-platform thresholds.

Ported from cpg-flow-rdrnaseq's scripts/check_multiqc.py, with three long-read differences:
- thresholds are looked up per SG by its sequencing platform ([qc_thresholds.<platform>.*]), not sequencing type
- mosdepth's metric names are mapped to the short-read names (MEDIAN_COVERAGE, PCT_20X), so one problem has one
  name in `cram_qc_flags` whichever technology it came from
- reads_mapped_percent is only taken from the whole-CRAM idxstats section: samtools stats runs on a sample of
  mapped reads, so its own reads_mapped_percent is always ~100%
"""

import json
from argparse import ArgumentParser, BooleanOptionalAction
from collections import defaultdict
from dataclasses import asdict
from datetime import datetime
from typing import Any

from loguru import logger

from rd_qc.utils import QC_DIRECTIONS, QcFlag, load_qc_thresholds, worst_qc_breach

from cpg_utils import config, to_path
from cpg_utils.config import dataset_for_access_level
from cpg_utils.slack import send_message
from metamist.graphql import gql, query

SG_PLATFORM_QUERY = gql(
    """
    query sgPlatforms($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups(technology: {eq: "long-read"}) {
                id
                platform
            }
        }
    }
    """
)

# mosdepth metric -> (short-read metric name, multiplier). MultiQC reports % >= 20x as a percentage; PCT_20X is a
# fraction, as Picard reports it for short read.
MOSDEPTH_METRIC_ALIASES: dict[str, tuple[str, float]] = {
    'median_coverage': ('MEDIAN_COVERAGE', 1.0),
    '20_x_pc': ('PCT_20X', 0.01),
}

# Metrics that are only valid from one MultiQC general stats section
METRIC_SECTIONS: dict[str, str] = {'reads_mapped_percent': 'lrs_idxstats'}


def get_sg_platforms(dataset: str) -> dict[str, str]:
    """{SG ID: sequencing platform} for the dataset's long-read SGs."""
    response = query(SG_PLATFORM_QUERY, variables={'dataset': dataset_for_access_level(dataset)})
    return {sg['id']: sg['platform'] for sg in response['project']['sequencingGroups']}


def general_stats_sections(multiqc_data: dict) -> dict[str, dict[str, dict[str, Any]]]:
    """
    MultiQC's general stats as {section name: {sample: {metric: value}}}, with mosdepth metrics renamed.

    Older MultiQC versions write the sections as an unnamed list, so those are named by position.
    """
    raw = multiqc_data['report_general_stats_data']
    sections = raw if isinstance(raw, dict) else {f'section_{i}': section for i, section in enumerate(raw)}

    for name, section in sections.items():
        if 'mosdepth' not in name:
            continue
        for val_by_metric in section.values():
            for metric, (alias, multiplier) in MOSDEPTH_METRIC_ALIASES.items():
                if metric in val_by_metric:
                    try:
                        val_by_metric[alias] = float(val_by_metric[metric]) * multiplier
                    except (TypeError, ValueError):
                        val_by_metric[alias] = val_by_metric[metric]
    return sections


def warn_unmatched_metrics(sections: dict[str, dict], platforms: set[str]) -> None:
    """Log a warning for any configured threshold metric MultiQC never surfaced, so a wrong key isn't silent."""
    present_metrics = {
        metric for section in sections.values() for val_by_metric in section.values() for metric in val_by_metric
    }
    for platform in sorted(platforms):
        configured = {metric for by_metric in load_qc_thresholds(platform).values() for metric in by_metric}
        if not configured:
            logger.warning(f'No qc_thresholds configured for platform={platform!r}; its SGs will not be checked.')
        for metric in sorted(configured - present_metrics):
            logger.warning(
                f'Configured threshold metric {metric!r} for platform={platform!r} not found in any MultiQC section; '
                'this threshold will not be checked.',
            )


def _check_sample(
    sample: str,
    section_name: str,
    val_by_metric: dict[str, Any],
    thresholds: dict[str, dict[str, dict[str, float]]],
    today: datetime,
) -> list[tuple[QcFlag, str]]:
    """The QC flags one sample's metrics in one section breach, each with its report line."""
    breaches: list[tuple[QcFlag, str]] = []
    for direction, metric_tiers in thresholds.items():
        sign = QC_DIRECTIONS[direction][0]
        for metric, tiers in metric_tiers.items():
            if metric not in val_by_metric:
                continue
            if metric in METRIC_SECTIONS and METRIC_SECTIONS[metric] not in section_name:
                continue
            try:
                val = float(val_by_metric[metric])
            except (TypeError, ValueError):
                logger.warning(
                    f'{sample}: metric {metric!r} has non-numeric value {val_by_metric[metric]!r}; '
                    'skipping threshold check.',
                )
                continue
            verdict = worst_qc_breach(val, tiers, direction)
            if verdict is None:
                logger.info(f'{sample}: {metric}={val:0.4f} within thresholds')
                continue
            severity, threshold = verdict
            icon = '❗' if severity == 'fail' else '⚠️'
            line = f'{icon} {metric}={val:0.4f}{sign}{threshold:0.4f} [{severity}]'
            logger.info(f'{sample}: {line}')
            flag = QcFlag(
                flag=metric,
                value=val,
                comparison=sign,
                threshold=threshold,
                section=section_name,
                date=today.isoformat(timespec='seconds'),
                ar_guid=config.try_get_ar_guid(),
                severity=severity,
            )
            breaches.append((flag, line))
    return breaches


def run(
    multiqc_json_path: str,
    html_url: str | None,
    dataset: str,
    title: str | None,
    send_to_slack: bool,
    output_json_path: str | None,
) -> dict[str, Any]:
    today = datetime.now()  # noqa: DTZ005

    with to_path(multiqc_json_path).open() as f:
        sections = general_stats_sections(json.load(f))
    logger.info(
        'report_general_stats_data: '
        + ', '.join(f'{name}={len(section)} samples' for name, section in sections.items())
    )

    platform_by_sg = get_sg_platforms(dataset)
    warn_unmatched_metrics(sections, set(platform_by_sg.values()))
    thresholds_by_platform = {platform: load_qc_thresholds(platform) for platform in set(platform_by_sg.values())}

    bad_lines_by_sample: dict[str, list[str]] = defaultdict(list)
    qc_flags_by_sample: dict[str, list[QcFlag]] = defaultdict(list)
    for section_name, section in sections.items():
        for sample, val_by_metric in section.items():
            sg_id = sample.split('|', 1)[0]
            platform = platform_by_sg.get(sg_id)
            if platform is None:
                logger.warning(f'{sample}: not a long-read SG in {dataset}; skipping threshold checks.')
                continue
            for flag, line in _check_sample(
                sample, section_name, val_by_metric, thresholds_by_platform[platform], today
            ):
                qc_flags_by_sample[sg_id].append(flag)
                bad_lines_by_sample[sample].append(line)

    report_title = title or 'Long-read MultiQC report'
    title_line = f'*[{dataset}]* <{html_url}|{report_title}>' if html_url else f'*[{dataset}]* {report_title}'
    messages = []
    if bad_lines_by_sample:
        n_fail = sum(1 for flags in qc_flags_by_sample.values() for f in flags if f.severity == 'fail')
        n_warn = sum(1 for flags in qc_flags_by_sample.values() for f in flags if f.severity == 'warn')
        messages.append(
            f'{title_line}. {len(bad_lines_by_sample)} samples flagged ({n_fail} failing, {n_warn} warnings):',
        )
        messages.extend(f'{sample}: ' + ', '.join(lines) for sample, lines in bad_lines_by_sample.items())
    else:
        messages.append(f'{title_line}. No samples flagged.')
    text = '\n'.join(messages)
    logger.info(text)

    result: dict[str, Any] = {
        'title': report_title,
        'dataset': dataset,
        'html_url': html_url,
        'sequencing_technology': 'long-read',
        'n_samples_flagged': len(qc_flags_by_sample),
        'qc_flags': {sg_id: [asdict(flag) for flag in flags] for sg_id, flags in qc_flags_by_sample.items()},
    }

    if output_json_path:
        with to_path(output_json_path).open('w') as f:
            json.dump(result, f, indent=2)

    if send_to_slack:
        try:
            send_message(text)
        except (OSError, ValueError):
            logger.exception('Failed to send Slack notification, continuing without it')

    return result


if __name__ == '__main__':
    parser = ArgumentParser(description='Check long-read MultiQC metrics against per-platform thresholds')
    parser.add_argument('--multiqc-json', required=True, help='Path to MultiQC JSON output')
    parser.add_argument('--html-url', help='MultiQC HTML URL')
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--title', help='Report title')
    parser.add_argument('--send-to-slack', action=BooleanOptionalAction, default=True, help='Post the summary')
    parser.add_argument('--output-json', help='Path to write structured QC flags JSON output')
    args = parser.parse_args()
    run(
        multiqc_json_path=args.multiqc_json,
        html_url=args.html_url,
        dataset=args.dataset,
        title=args.title,
        send_to_slack=args.send_to_slack,
        output_json_path=args.output_json,
    )
