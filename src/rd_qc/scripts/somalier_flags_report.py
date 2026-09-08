"""
Queries Metamist for all Somalier flags across a dataset's sequencing groups
and renders them into a report using the somalier_flags_report.html.jinja template.
"""

import re
from argparse import ArgumentParser
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from time import perf_counter

import jinja2
from loguru import logger

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, dataset_for_access_level
from cpg_utils.metamist_registration import create_new
from cpg_utils.slack import send_message
from metamist.graphql import gql, query

from rd_qc.utils import SomalierFlag

STAGE_NAME = 'GenerateSomalierFlagsReport'
JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'

DATASET_SGS_QUERY = gql(
    """
    query datasetSgs($dataset: String!, $seqType: String!, $seqTech: String!) {
        project(name: $dataset) {
            sequencingGroups(type: {eq: $seqType}, technology: {eq: $seqTech}) {
                id
                meta
                type
            }
        }
    }
    """
)

SGS_INFO_QUERY = gql(
    """
    query sgInfo($sgIds: [String!]!) {
        sequencingGroups(id: {in_: $sgIds}) {
            id
            meta
            type
            technology
            platform
            assays {
                id
                meta
            }
            sample {
                id
                externalIds
                type
                participant {
                    id
                    externalIds
                    families {
                        id
                        externalIds
                    }
                }
            }
        }
    }
    """
)

EXISTING_ANALYSES_QUERY = gql(
    """
    query existingAnalyses($dataset: String!, $metaFilter: JSON!) {
        project(name: $dataset) {
            analyses(type: {eq: "web"}, meta: $metaFilter) {
                id
                outputs
                timestampCompleted
                meta
            }
        }
    }
    """
)

@dataclass
class SGInfo:
    sg_id: str
    sg_type: str
    sg_technology: str
    sg_platform: str
    crams: list[str]
    fastq_pairs: list[tuple[str, str]]
    other_reads: list[str]
    sample_external_id: str
    sample_type: str
    participant_external_id: str
    family_external_id: str


@dataclass(frozen=True)
class SGReport:
    """All QC flags for a single sequencing group, plus its metadata."""

    sg_info: SGInfo
    somalier_flags: list[SomalierFlag]

def _has_active(flags: list[SomalierFlag]) -> bool:
    return any(not f.resolved for f in flags)

# ---------------------------------------------------------------------------
# Metamist querying
# ---------------------------------------------------------------------------
def _primary_external_id(obj: dict | None) -> str:
    """Metamist keys the primary external ID under ''. Fall back gracefully."""
    ext = (obj or {}).get('externalIds') or {}
    return ext.get('') or next(iter(ext.values()), '')


def _basename(entry: dict | str | None) -> str | None:
    """Pull a file basename from a reads entry (dict or path string)."""
    if isinstance(entry, dict):
        return entry.get('basename') or ((entry.get('location') or '').rsplit('/', 1)[-1] or None)
    if isinstance(entry, str):
        return entry.rsplit('/', 1)[-1] or None
    return None


def _extract_reads(assays: list[dict]) -> tuple[list[str], list[tuple[str, str]], list[str]]:
    """Group an SG's assay read files into (crams, fastq_pairs, other).

    Each assay's ``meta.reads`` is a list of file entries; ``meta.reads_type``
    tells us whether they're fastq (R1/R2 pairs) or aligned (bam/cram).
    """
    crams: list[str] = []
    fastq_pairs: list[tuple[str, str]] = []
    other: list[str] = []

    for assay in assays:
        meta = assay.get('meta') or {}
        reads = meta.get('reads')
        reads_type = (meta.get('reads_type') or '').lower()
        entries = reads if isinstance(reads, list) else ([reads] if reads else [])
        names = [n for n in (_basename(e) for e in entries) if n]
        if not names:
            continue

        if reads_type == 'fastq':
            # A fastq assay is one (or more) R1/R2 pair(s); pair sequentially.
            for i in range(0, len(names), 2):
                r2 = names[i + 1] if i + 1 < len(names) else ''
                fastq_pairs.append((names[i], r2))
        elif reads_type in ('bam', 'cram'):
            crams.extend(names)
        else:
            # Unknown reads_type — classify by extension.
            for name in names:
                low = name.lower()
                if low.endswith(('.cram', '.bam')):
                    crams.append(name)
                elif low.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')):
                    fastq_pairs.append((name, ''))
                else:
                    other.append(name)

    return crams, fastq_pairs, other


def get_sg_infos(sg_ids: list[str]) -> dict[str, SGInfo]:
    """Query Metamist for detailed SG info, keyed by SG id."""
    if not sg_ids:
        logger.warning('No SGs flagged')
        return {}
    logger.info(f'Querying Metamist for detailed info on {len(sg_ids)} SG(s): {", ".join(sg_ids)}')
    started = perf_counter()
    response = query(SGS_INFO_QUERY, variables={'sgIds': sg_ids})
    logger.info(f'Received SG info for {len(sg_ids)} SG(s) in {perf_counter() - started:.1f}s')
    infos: dict[str, SGInfo] = {}
    for sg in response['sequencingGroups']:
        sample = sg.get('sample') or {}
        participant = sample.get('participant') or {}
        families = participant.get('families') or []
        family = families[0] if families else None

        crams, fastq_pairs, other = _extract_reads(sg.get('assays', []))

        infos[sg['id']] = SGInfo(
            sg_id=sg['id'],
            sg_type=sg.get('type') or '',
            sg_technology=sg.get('technology') or '',
            sg_platform=sg.get('platform') or '',
            crams=crams,
            fastq_pairs=fastq_pairs,
            other_reads=other,
            sample_external_id=_primary_external_id(sample),
            sample_type=sample.get('type') or '',
            participant_external_id=_primary_external_id(participant),
            family_external_id=_primary_external_id(family),
        )
    return infos


def get_previous_analysis(dataset: str, meta_filter: dict) -> dict | None:
    """Query Metamist for existing web analyses matching a meta filter and take the most recent one, if any."""
    response = query(EXISTING_ANALYSES_QUERY, variables={'dataset': dataset, 'metaFilter': meta_filter})
    existing_analyses = response['project']['analyses']
    if not existing_analyses:
        return None
    existing_analyses.sort(key=lambda a: a.get('timestampCompleted') or '', reverse=True)
    if len(existing_analyses) == 1:
        # Only one analysis is the current one, so there are no previous analyses to compare against.
        return None
    previous_analysis = existing_analyses[1]  # The second most recent analysis is the previous one
    logger.info(f'Found previous analysis {previous_analysis["id"]} from {previous_analysis["timestampCompleted"]}')
    if not previous_analysis['meta'].get('summary'):
        logger.warning(f'Previous analysis {previous_analysis["id"]} has no summary in meta; skipping')
        return None
    return previous_analysis


def collect_somalier_flags(sequencing_groups: list[dict]) -> list[dict]:
    """Extract Somalier flags from each sequencing group's metadata."""
    results = []
    for sg in sequencing_groups:
        meta = sg.get('meta') or {}
        results.append(
            {
                'id': sg['id'],
                'somalier_flags': [SomalierFlag(**flag) for flag in meta.get('somalier_flags', [])]
            }
        )
    return results

def summarise_flags(sg_data: list[dict]) -> dict:
    """Dataset-wide, flag-centric summary counts for the header cards."""
    all_flags: list[SomalierFlag] = [f for sg in sg_data for f in sg['somalier_flags']]

    active = [f for f in all_flags if not f.resolved]
    active_sex_inference_flags = sum(1 for f in active if f.category == 'sex_inference_mismatch')
    active_self_relatedness_flags = sum(1 for f in active if f.category == 'self_relatedness_mismatch')
    active_relatedness_flags = sum(1 for f in all_flags if f.category == 'relatedness_mismatch')
    sgs_affected = sum(1 for sg in sg_data if _has_active(sg['somalier_flags']))

    return {
        'total_sgs': len(sg_data),
        'active_flags': len(active),
        'active_sex_inference_flags': active_sex_inference_flags,
        'active_self_relatedness_flags': active_self_relatedness_flags,
        'active_relatedness_flags': active_relatedness_flags,
        'sgs_affected': sgs_affected,
        'resolved_flags': sum(1 for f in all_flags if f.resolved),
    }

def render_report(dataset: str, reports: list[SGReport], summary: dict) -> str:
    """Render the Somalier flags report HTML using the Jinja template."""
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
        autoescape=jinja2.select_autoescape(['html', 'xml']),
    )
    template = env.get_template('somalier_flags_report.html.jinja')
    return template.render(dataset=dataset, reports=reports, summary=summary)


def construct_summary_message(  # noqa: PLR0917
    dataset: str,
    out_html_url: str,
    somalier_url: str,
    seq_type: str,
    seq_tech: str,
    summary: dict,
    previous_analysis: dict | None,
):
    """Construct a Slack message with a concise summary and a link to the report."""
    pass


def main(
    dataset: str,
    output: str,
    timestamped_output: str,
    out_html_url: str,
    somalier_url: str,
) -> None:
    """Query Metamist for Somalier flags and generate a Somalier flags HTML report."""

    dataset = dataset_for_access_level(dataset)
    seq_type = config_retrieve(['workflow', 'sequencing_type'])
    seq_tech = config_retrieve(['workflow', 'sequencing_technology'])

    logging_prefix = f'{dataset} ({seq_type} | {seq_tech})'

    logger.info(f'{logging_prefix} :: Querying Metamist for Somalier flags')
    started = perf_counter()
    response = query(DATASET_SGS_QUERY, variables={'dataset': dataset, 'seqType': seq_type, 'seqTech': seq_tech})
    sequencing_groups = response['project']['sequencingGroups']
    logger.info(
        f'{logging_prefix} :: Found {len(sequencing_groups)} sequencing groups in {perf_counter() - started:.1f}s'
    )

    sg_data = collect_somalier_flags(sequencing_groups)
    summary = summarise_flags(sg_data)

    # Only fetch rich SG metadata for SGs that actually have flags to report.
    flagged = [sg for sg in sg_data if sg['somalier_flags']]
    logger.info(f'{logging_prefix} :: {len(flagged)} sequencing groups have Somalier flags.')

    infos = get_sg_infos([sg['id'] for sg in flagged])
    reports = [
        SGReport(
            sg_info=infos[sg['id']],
            somalier_flags=sg['somalier_flags'],
        )
        for sg in flagged
        if sg['id'] in infos
    ]

    logger.info(f'{logging_prefix} :: Rendering report for {len(reports)} flagged SG(s)')
    started = perf_counter()
    html = render_report(dataset, reports, summary=summary)
    logger.info(f'{logging_prefix} :: Rendered report in {perf_counter() - started:.1f}s')

    with to_path(output).open('w') as f:
        f.write(html)
    logger.info(f'{logging_prefix} :: Wrote SG QC report to {output}')

    with to_path(timestamped_output).open('w') as f:
        f.write(html)
    logger.info(f'{logging_prefix} :: Wrote timestamped SG QC report to {timestamped_output}')

    # Register results in Metamist manually to capture all dataset SGs in scope, not just the input_cohorts SGs
    meta = {
        'stage': STAGE_NAME,
        'dataset': dataset,
        'sequencing_type': seq_type,
        'sequencing_technology': seq_tech,
        'summary': summary,
    }

    create_new(
        project=dataset,
        output=timestamped_output,
        analysis_type='web',
        sgs=[sg['id'] for sg in sequencing_groups],
        meta=meta,
    )
    logger.info(f'{logging_prefix} :: Registered web analysis for {len(sequencing_groups)} SG(s)')

    meta.pop('summary')
    construct_summary_message(
        dataset,
        out_html_url,
        somalier_url,
        seq_type,
        seq_tech,
        summary,
        get_previous_analysis(dataset, meta),
    )



if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Metamist dataset/project name')
    parser.add_argument('--fixed-output', required=True, help='Path to write the HTML report')
    parser.add_argument('--timestamped-output', required=True, help='Path to write the timestamped HTML report')
    parser.add_argument('--html-url', required=True, help='Clickable URL for the Somalier flags HTML report')
    parser.add_argument('--somalier-url', required=True, help='Clickable URL for the original Somalier report')
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        output=args.fixed_output,
        timestamped_output=args.timestamped_output,
        out_html_url=args.html_url,
        somalier_url=args.somalier_url
    )
