"""
Batch jobs to run MultiQC over long-read CRAM QC, check thresholds, and record QC flags.

Ported from cpg-flow-rdrnaseq's jobs/multiqc.py: inputs are the dataset's registered `qc` analyses rather than
stage inputs, so SGs whose QC ran in an earlier workflow run are still included.
"""

from collections import defaultdict
from os.path import basename

from hailtop.batch import Batch, Resource, ResourceFile
from hailtop.batch.job import Job
from loguru import logger

from cpg_flow import resources, targets, utils
from cpg_utils import Path, config, hail_batch
from metamist.graphql import gql, query

QC_ANALYSES_QUERY = gql(
    """
    query QcAnalyses($dataset: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: "qc"}, status: {eq: COMPLETED}) {
                output
                meta
                timestampCompleted
                sequencingGroups {
                    id
                }
            }
        }
    }
""",
)

# The per-SG long-read QC stages. Each registers one `qc` analysis with meta.stage set to its name.
LRS_QC_STAGES = {'LongReadQcSamtoolsStats', 'LongReadQcVerifyBamId'}

# mosdepth is not run here: its outputs come with the long-read data and are ingested as `qc` analyses with
# meta.tool set to this.
MOSDEPTH_TOOL = 'mosdepth'

# Files written next to a registered output, found by name: {registered suffix: [companion suffixes]}
COMPANION_SUFFIXES = {
    '.samtools-stats': ['.lrs_mapped_mqc.tsv'],
    '.mosdepth.global.dist.txt': ['.mosdepth.summary.txt'],
    '.mosdepth.summary.txt': ['.mosdepth.global.dist.txt'],
}


def _qc_source(meta: dict) -> str | None:
    """Which long-read QC source a `qc` analysis comes from, or None if it isn't one MultiQC reads."""
    if meta.get('stage') in LRS_QC_STAGES:
        return meta['stage']
    if meta.get('tool') == MOSDEPTH_TOOL:
        return MOSDEPTH_TOOL
    return None


def _get_registered_qc_paths(dataset_name: str) -> dict[str, list[str]]:
    """
    Query Metamist for the dataset's registered long-read QC outputs: those of the long-read QC stages, and the
    ingested mosdepth outputs.

    Returns {SG ID: [output paths]}, taking the latest analysis per SG per source per file type, so a source
    registering several files (e.g. mosdepth's dist and summary) keeps all of them.
    """
    result = query(QC_ANALYSES_QUERY, variables={'dataset': dataset_name})

    latest: dict[tuple[str, str, str], tuple[str, str]] = {}
    for analysis in result.get('project', {}).get('analyses', []):
        source = _qc_source(analysis.get('meta') or {})
        if source is None:
            continue
        sg_ids = [sg['id'] for sg in analysis.get('sequencingGroups', [])]
        if len(sg_ids) != 1:
            continue
        output = analysis['output']
        file_type = basename(output).split('.', 1)[-1]
        key = (sg_ids[0], source, file_type)
        ts = analysis.get('timestampCompleted', '')
        if key not in latest or ts > latest[key][0]:
            latest[key] = (ts, output)

    sg_paths: dict[str, list[str]] = defaultdict(list)
    for (sg_id, _source, _file_type), (_ts, output) in latest.items():
        sg_paths[sg_id].append(output)
    return dict(sg_paths)


def _add_companion_paths(sg_paths: dict[str, list[str]]) -> None:
    """Add the files written next to each registered output (see COMPANION_SUFFIXES), if they exist."""
    for paths in sg_paths.values():
        for path in list(paths):
            for suffix, companions in COMPANION_SUFFIXES.items():
                if not path.endswith(suffix):
                    continue
                for companion in companions:
                    companion_path = path.removesuffix(suffix) + companion
                    if companion_path not in paths and utils.can_reuse(companion_path):
                        paths.append(companion_path)


def lrs_multiqc(
    dataset: targets.Dataset,
    outputs: dict[str, Path],
    out_html_url: str,
    label: str,
    job_attrs: dict,
) -> list[Job]:
    """
    Run MultiQC, then check thresholds, then record flags in Metamist.

    All three jobs are returned, so a stage depending on this one waits for the flags to be written.
    `label` sets the SG meta key the flags go into: 'CRAM' -> 'cram_qc_flags'.
    """
    batch_instance = hail_batch.get_batch()
    tmp_prefix = dataset.tmp_prefix() / 'multiqc' / 'long_read'
    dataset_name = config.dataset_for_access_level(dataset.name)
    sequencing_group_id_map = dataset.rich_id_map()

    sg_paths = _get_registered_qc_paths(dataset_name)
    _add_companion_paths(sg_paths)
    logger.info(f'Long-read MultiQC for {dataset.name}: {len(sg_paths)} SGs with registered QC outputs')

    all_paths = {p for sg_files in sg_paths.values() for p in sg_files}
    if not all_paths:
        logger.warning(f'No long-read QC outputs found for dataset {dataset.name}, skipping MultiQC')
        return []

    title = f'MultiQC [Long-read {label}]'

    mqc_j = batch_instance.new_bash_job(title, job_attrs | {'tool': 'MultiQC'})
    mqc_j.image(config.config_retrieve(['images', 'multiqc']))
    mqc_j.cpu(8)
    mqc_j.storage('20Gi')

    sg_id_mapping_file_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_rename-sg-map.tsv'
    if not config.config_retrieve(['workflow', 'dry_run'], False):
        with sg_id_mapping_file_path.open('w') as fh:
            for sgid, new_sgid in sequencing_group_id_map.items():
                fh.write('\t'.join([sgid, new_sgid]) + '\n')
    sg_id_mapping_file = batch_instance.read_input(sg_id_mapping_file_path)

    qc_inputs = batch_instance.read_input_group(**{basename(p): p for p in sorted(all_paths)})

    endings = ', '.join(['.mosdepth', '.samtools-stats', '.verify-bamid.selfSM', '.lrs_mapped_mqc'])

    mqc_j.command(
        f"""
        multiqc -f $(dirname {qc_inputs}) -o output \\
        --replace-names {sg_id_mapping_file} \\
        --title "{title} for dataset <b>{dataset.name}</b>" \\
        --filename report.html \\
        --cl-config "extra_fn_clean_exts: [{endings}]" \\
        --cl-config "max_table_rows: 10000"

        cp output/report.html {mqc_j.html}
        cp output/report_data/multiqc_data.json {mqc_j.json}
        """
    )
    mqc_j.command(f'echo "HTML URL: {out_html_url}"')

    batch_instance.write_output(mqc_j.html, str(outputs['html']))
    batch_instance.write_output(mqc_j.html, str(outputs['latest']))
    batch_instance.write_output(mqc_j.json, str(outputs['json']))
    all_jobs: list[Job] = [mqc_j]

    check_j = _check_report_job(
        b=batch_instance,
        multiqc_json_file=mqc_j.json,
        multiqc_html_url=out_html_url,
        rich_id_map=sequencing_group_id_map,
        dataset_name=dataset.name,
        label=label,
        out_checks_path=outputs['checks'],
        job_attrs=job_attrs,
    )
    check_j.depends_on(mqc_j)
    all_jobs.append(check_j)

    record_j = _record_qc_flags_job(
        b=batch_instance,
        dataset_name=dataset.name,
        label=label,
        sg_id_mapping_file=sg_id_mapping_file,
        check_multiqc_json_file=check_j.output,
        job_attrs=job_attrs,
    )
    record_j.depends_on(check_j)
    all_jobs.append(record_j)

    return all_jobs


def _rich_sequencing_group_id_seds(
    rich_id_map: dict[str, str],
    file_names: list[str | Resource | ResourceFile],
) -> str:
    """Extend sequencing group IDs in files with external IDs via sed."""
    cmd = ''
    for sgid, rich_sgid in rich_id_map.items():
        for fname in file_names:
            cmd += f"sed -iBAK 's/{sgid}/{rich_sgid}/g' {fname}\n"
    return cmd


def _check_report_job(
    b: Batch,
    multiqc_json_file: Resource | ResourceFile,
    dataset_name: str,
    multiqc_html_url: str,
    label: str,
    out_checks_path: Path,
    rich_id_map: dict[str, str] | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """Check the MultiQC JSON against the platform thresholds, and send a Slack summary if configured."""
    title = f'MultiQC [Long-read {label}]'
    check_j = b.new_job(f'{title} check', (job_attrs or {}) | {'tool': 'python'})
    resources.STANDARD.set_resources(j=check_j, ncpu=2)
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))

    send_to_slack = config.config_retrieve(['lrs_qc', 'multiqc', 'send_to_slack'], True)

    cmd = f"""\
    {_rich_sequencing_group_id_seds(rich_id_map, [multiqc_json_file]) if rich_id_map else ''}

    python3 -m rd_qc.scripts.check_lrs_qc \\
    --multiqc-json {multiqc_json_file} \\
    --html-url {multiqc_html_url} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --{'' if send_to_slack else 'no-'}send-to-slack \\
    --output-json {check_j.output}
    """

    check_j.command(cmd)
    b.write_output(check_j.output, str(out_checks_path))
    return check_j


def _record_qc_flags_job(
    b: Batch,
    dataset_name: str,
    label: str,
    sg_id_mapping_file: Resource | ResourceFile,
    check_multiqc_json_file: Resource | ResourceFile,
    job_attrs: dict | None = None,
) -> Job:
    """Record QC flags in Metamist, under the `{label.lower()}_qc_flags` SG meta key."""
    record_j = b.new_job('Record long-read QC flags', (job_attrs or {}) | {'tool': 'python'})
    resources.STANDARD.set_resources(j=record_j, ncpu=2)
    record_j.image(config.config_retrieve(['workflow', 'driver_image']))

    cmd = f"""\
    python3 -m rd_qc.scripts.record_lrs_qc_flags \\
    --dataset {dataset_name} \\
    --label {label} \\
    --qc-flags-json {check_multiqc_json_file} \\
    --sequencing-group-ids-map {sg_id_mapping_file}
    """

    record_j.command(cmd)
    return record_j
