"""
This script parses "somalier relate" (https://github.com/brentp/somalier) outputs,
and returns a report whether sex and pedigree matches the provided PED file.

Script can send a report to a Slack channel. To enable that, set SLACK_TOKEN
and SLACK_CHANNEL environment variables, and add "Seqr Loader" app into
a channel with:

/invite @Seqr Loader
"""

import contextlib
import json
from argparse import ArgumentParser
from dataclasses import asdict
from typing import Any

import pandas as pd
from loguru import logger
from peddy import Ped, Sample

from rd_qc.utils import SomalierRelatednessFlag, SomalierSexInferenceFlag

from cpg_flow.metamist import get_metamist
from cpg_utils import config, slack, to_path
from cpg_utils.metamist_registration import create_new

_messages: list[str] = []


def info(msg):
    _messages.append(msg)
    logger.info(msg)


def warning(msg):
    _messages.append(msg)
    logger.warning(msg)


def error(msg):
    _messages.append(msg)
    logger.error(msg)


def _format_mismatch_line(s1, s2, expected_ped_s1, expected_ped_s2, expected_rel, inferred_rel, row) -> str:
    fam1 = expected_ped_s1.family_id if expected_ped_s1 else None
    fam2 = expected_ped_s2.family_id if expected_ped_s2 else None
    if fam1 == fam2:
        line = f'{fam1}: {s1} - {s2}'
    else:
        line = s1 + (f' ({fam1})' if fam1 and fam1 != s1 else '')
        line += ' - '
        line += s2 + (f' ({fam2})' if fam2 and fam2 != s2 else '')
    return (
        f'{line}, '
        f'provided: "{expected_rel}", '
        f'inferred: "{inferred_rel}", '
        f'kin={row["relatedness"]}, '
        f'ibs0={row["ibs0"]}, '
        f'ibs2={row["ibs2"]}'
    )


def _check_sex(samples_df) -> dict[str, SomalierSexInferenceFlag]:
    info('*Inferred vs. reported sex:*')
    samples_df.sex = samples_df.sex.apply(lambda x: {1: 'male', 2: 'female'}.get(x, 'unknown'))
    samples_df.original_pedigree_sex = samples_df.original_pedigree_sex.apply(lambda x: {'-9': 'unknown'}.get(x, x))
    samples_df['x_het_ratio'] = samples_df.X_het / samples_df.X_hom_alt.replace(0, 1)  # avoid division by zero
    samples_df['x_depth_ratio'] = 2 * samples_df.X_depth_mean / samples_df.gt_depth_mean
    samples_df['y_depth_ratio'] = 2 * samples_df.Y_depth_mean / samples_df.gt_depth_mean
    missing_inferred_sex = samples_df.sex == 'unknown'
    missing_provided_sex = samples_df.original_pedigree_sex == 'unknown'
    mismatching_female = (samples_df.sex == 'female') & (samples_df.original_pedigree_sex == 'male')
    mismatching_male = (samples_df.sex == 'male') & (samples_df.original_pedigree_sex == 'female')
    mismatching_sex = mismatching_female | mismatching_male
    mismatching_other = (
        (samples_df.sex != samples_df.original_pedigree_sex) & (~mismatching_female) & (~mismatching_male)
    )
    matching_sex = ~mismatching_sex & ~mismatching_other

    sex_mismatches_by_sgid: dict[str, SomalierSexInferenceFlag] = {
        row.sample_id: SomalierSexInferenceFlag(
            category='sex_inference_mismatch',
            provided=str(row.original_pedigree_sex),
            inferred=str(row.sex),
            mean_depth=float(row.gt_depth_mean),
            x_het_ratio=float(row.x_het_ratio),
            x_depth_ratio=float(row.x_depth_ratio),
            y_depth_ratio=float(row.y_depth_ratio),
            x_sites=int(row.X_n),
            p_middling_ab=float(row.p_middling_ab),
        )
        for _, row in samples_df[mismatching_sex].iterrows()
    }

    def _print_stats(df_filter) -> None:
        for _, row_ in samples_df[df_filter].iterrows():
            info(
                f' {row_.sample_id} ('
                f'provided: {row_.original_pedigree_sex}, '
                f'inferred: {row_.sex}, '
                f'X het ratio: {row_.x_het_ratio:.2f})'
            )

    if mismatching_sex.any():
        info(f'❗ {len(samples_df[mismatching_sex])}/{len(samples_df)} PED samples with mismatching sex:')
        _print_stats(mismatching_sex)
    if missing_provided_sex.any():
        info(f'⚠️ {len(samples_df[missing_provided_sex])}/{len(samples_df)} samples with missing provided sex:')
        _print_stats(missing_provided_sex)
    if missing_inferred_sex.any():
        info(f'⚠️ {len(samples_df[missing_inferred_sex])}/{len(samples_df)} samples with failed inferred sex:')
        _print_stats(missing_inferred_sex)
    inferred_cnt = len(samples_df[~missing_inferred_sex])
    matching_cnt = len(samples_df[matching_sex])
    info(
        f'✅ Sex inferred for {inferred_cnt}/{len(samples_df)} samples, matching '
        f'for {matching_cnt if matching_cnt != inferred_cnt else "all"} samples.',
    )
    info('')

    return sex_mismatches_by_sgid


def _report_relatedness_findings(
    unrelated_to_related: list[str],
    related_to_unrelated: list[str],
) -> None:
    if unrelated_to_related:
        info(
            f'⚠️ Found {len(unrelated_to_related)} '
            f'sample pair(s) that are provided as unrelated, are inferred as '
            f'related:',
        )
        for i, pair in enumerate(unrelated_to_related):
            info(f' {i + 1}. {pair}')
    if related_to_unrelated:
        info(
            f'❗ Found {len(related_to_unrelated)} sample pair(s) '
            f'that are provided as related, but inferred as unrelated:',
        )
        for i, pair in enumerate(related_to_unrelated):
            info(f' {i + 1}. {pair}')
    if not unrelated_to_related and not related_to_unrelated:
        info('✅ Inferred pedigree matches for all provided related pairs.')
    info('')


def _check_relatedness(
    pairs_df,
    expected_ped: Ped,
    inferred_ped: Ped,
    bad_ids: list,
) -> dict[str, list[SomalierRelatednessFlag]]:
    info('*Relatedness:*')
    expected_ped_sample_by_id: dict[str, Sample] = {s.sample_id: s for s in expected_ped.samples()}
    inferred_ped_sample_by_id: dict[str, Sample] = {s.sample_id: s for s in inferred_ped.samples()}

    mismatching_unrelated_to_related = []
    mismatching_related_to_unrelated = []

    relatedness_flags_by_sg_id: dict[str, list[SomalierRelatednessFlag]] = {}
    for idx, row in pairs_df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        if s1 in bad_ids or s2 in bad_ids:
            continue

        expected_ped_s1 = expected_ped_sample_by_id.get(s1, {})
        expected_ped_s2 = expected_ped_sample_by_id.get(s2, {})
        inferred_ped_s1 = inferred_ped_sample_by_id.get(s1, {})
        inferred_ped_s2 = inferred_ped_sample_by_id.get(s2, {})
        with contextlib.redirect_stderr(None), contextlib.redirect_stdout(None):
            if expected_ped_s1 and expected_ped_s2:
                expected_rel = expected_ped.relation(expected_ped_s1, expected_ped_s2)
            else:
                expected_rel = 'unknown'
            if inferred_ped_s1 and inferred_ped_s2:
                inferred_rel = inferred_ped.relation(inferred_ped_s1, inferred_ped_s2)
            else:
                inferred_rel = 'unknown'

        if inferred_rel != expected_rel:
            line = _format_mismatch_line(s1, s2, expected_ped_s1, expected_ped_s2, expected_rel, inferred_rel, row)
            # peddy .samples() yields Sample objects (attribute access), but the
            # dict lookup above falls back to {} when a sample is missing, so guard
            # both cases with getattr.
            family_external_id = (
                getattr(expected_ped_s1, 'family_id', None) or getattr(expected_ped_s2, 'family_id', None) or 'unknown'
            )
            if s1 not in relatedness_flags_by_sg_id:
                relatedness_flags_by_sg_id[s1] = []
            relatedness_flags_by_sg_id[s1].append(
                SomalierRelatednessFlag(
                    category='relatedness_mismatch',
                    sg_id_1=s1,
                    sg_id_2=s2,
                    family_external_id=family_external_id,
                    expected_relationship=expected_rel,
                    inferred_relationship=inferred_rel,
                    relatedness=row['relatedness'],
                    ibs0=row['ibs0'],
                    ibs2=row['ibs2'],
                ),
            )

            if (expected_rel == 'unknown' and inferred_rel != 'unknown') or (
                expected_rel == 'unrelated' and inferred_rel != 'unrelated'
            ):
                if row['relatedness'] > 0.1:  # noqa: PLR2004
                    mismatching_unrelated_to_related.append(line)
            else:
                mismatching_related_to_unrelated.append(line)

        pairs_df.loc[idx, 'provided_rel'] = expected_rel
        pairs_df.loc[idx, 'inferred_rel'] = inferred_rel

    _report_relatedness_findings(mismatching_unrelated_to_related, mismatching_related_to_unrelated)

    return relatedness_flags_by_sg_id


def run(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    expected_ped_fpath: str,
    title: str,
    sg_ids: list[str],
    output_pairs: str,
    output_samples: str,
    output_html: str,
    html_url: str,
    base_output_html: str,
    dataset: str,
    output_json: str | None = None,
):
    """Report pedigree inconsistencies, given somalier outputs."""

    dataset = get_metamist().get_metamist_proj(dataset)

    logger.info(somalier_samples_fpath)
    samples_df = pd.read_csv(somalier_samples_fpath, delimiter='\t')
    pairs_df = pd.read_csv(somalier_pairs_fpath, delimiter='\t')
    with to_path(somalier_samples_fpath).open() as f:
        inferred_ped = Ped(f)
    with to_path(expected_ped_fpath).open() as f:
        expected_ped = Ped(f)

    bad = samples_df.gt_depth_mean == 0.0
    if bad.any():
        warning(
            f'⚠️ Excluded {len(samples_df[bad])}/{len(samples_df)} samples with zero '
            f'mean GT depth from pedigree/sex checks: {", ".join(samples_df[bad].sample_id)}',
        )
        info('')
    bad_ids = list(samples_df[bad].sample_id)  # for checking in pairs_df
    samples_df = samples_df[~bad]

    sex_mismatches_by_sgid = _check_sex(samples_df)
    relatedness_flags_by_sg_id = _check_relatedness(
        pairs_df,
        expected_ped,
        inferred_ped,
        bad_ids,
    )

    print_contents(
        samples_df,
        pairs_df,
        somalier_samples_fpath,
        somalier_pairs_fpath,
    )

    if dataset and html_url:
        title = f'*[{dataset}]* <{html_url}|{title or "Somalier pedigree report"}>'
    elif not title:
        title = 'Somalier pedigree report'
    text = '\n'.join([title, *_messages])

    if config.config_retrieve(['somalier_pedigree', 'send_to_slack'], default=True):
        slack.send_message(text)

    all_flags_by_sg_id: dict[str, list[SomalierSexInferenceFlag | SomalierRelatednessFlag]] = {}
    for sg_id, flag in sex_mismatches_by_sgid.items():
        all_flags_by_sg_id[sg_id] = [flag]
    for sg_id, flags in relatedness_flags_by_sg_id.items():
        if sg_id not in all_flags_by_sg_id:
            all_flags_by_sg_id[sg_id] = []
        all_flags_by_sg_id[sg_id].extend(flags)

    result: dict[str, Any] = {
        'dataset': dataset,
        'html_url': html_url,
        'n_samples_flagged': len([flags for flags in all_flags_by_sg_id.values() if flags]),
        'relatedness_flags': {sg_id: [asdict(flag) for flag in flags] for sg_id, flags in all_flags_by_sg_id.items()},
    }
    if output_json:
        with to_path(output_json).open('w') as f:
            json.dump(result, f, indent=2)

    with to_path(output_samples).open('w') as f:
        f.write(to_path(somalier_samples_fpath).read_text())
    with to_path(output_pairs).open('w') as f:
        f.write(to_path(somalier_pairs_fpath).read_text())

    # Now create the web analysis for the whole dataset
    create_new(
        project=dataset,
        output=output_html,
        analysis_type='web',
        sgs=sg_ids,
        meta={'stage': 'SomalierPedigreeCheck'},
        secondary={
            'base_html': base_output_html,
            'samples': output_samples,
            'pairs': output_pairs,
            'json': output_json,
        },
    )
    logger.info(f'Registered web analysis for {dataset} at {output_html}')

    return result


def print_contents(
    samples_df,
    pairs_df,
    somalier_samples_fpath,
    somalier_pairs_fpath,
):
    if len(samples_df) < 400:  # noqa: PLR2004
        samples_str = samples_df.to_string()
        logger.info(f'Somalier results, samples (based on {somalier_samples_fpath}):\n{samples_str}\n')
    if len(pairs_df) < 400:  # noqa: PLR2004
        pairs_str = pairs_df[
            [
                '#sample_a',
                'sample_b',
                'relatedness',
                'ibs0',
                'ibs2',
                'n',
                'expected_relatedness',
            ]
        ].to_string()
        logger.info(f'Somalier results, sample pairs (based on {somalier_pairs_fpath}):\n{pairs_str}\n')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '--somalier-samples',
        required=True,
        help='Path to somalier {prefix}.samples.tsv output file',
    )
    parser.add_argument(
        '--somalier-pairs',
        required=True,
        help='Path to somalier {prefix}.pairs.tsv output file',
    )
    parser.add_argument(
        '--ped',
        required=True,
        help='Path to PED file with expected pedigree',
    )
    parser.add_argument('--title', required=True, help='Report title')
    parser.add_argument('--html-url', help='Web-accessible HTML path (namespaced by AR GUID)')
    parser.add_argument('--dataset', help='Dataset name')
    parser.add_argument('--sg-ids', nargs='+', required=True, help='space-separated SG IDs')
    parser.add_argument('--output-pairs', required=True)
    parser.add_argument('--output-samples', required=True)
    parser.add_argument('--output-html', required=True, help='gs:// path to HTML (namespaced by AR GUID)')
    parser.add_argument('--base-output-html', required=True, help='gs:// path to HTML (fixed, not namespaced)')
    parser.add_argument('--output-json', required=False, help='JSON output file for results')
    args = parser.parse_args()
    run(
        somalier_samples_fpath=args.somalier_samples,
        somalier_pairs_fpath=args.somalier_pairs,
        expected_ped_fpath=args.ped,
        html_url=args.html_url,
        dataset=args.dataset,
        title=args.title,
        sg_ids=args.sg_ids,
        output_pairs=args.output_pairs,
        output_samples=args.output_samples,
        output_html=args.output_html,
        base_output_html=args.base_output_html,
        output_json=args.output_json,
    )
