"""
Jobs for somalier relate — used by both identity checks and pedigree checks.
"""

from hailtop.batch.job import BashJob

from cpg_utils import Path, config, hail_batch


def _get_out_html_url(dataset_name: str, html_path: Path | str) -> str:
    """
    Convert a gs:// web-bucket path to the http(s) web URL.
    """
    # Important - strip -test from dataset suffix before constructing the web URL
    dataset_name = dataset_name.removesuffix('-test')
    return str(html_path).replace(
        config.config_retrieve(['storage', dataset_name, 'web']),
        config.config_retrieve(['storage', dataset_name, 'web_url']),
    )


def self_relatedness_jobs(
    dataset_name: str,
    participant_id: str,
    somalier_paths: dict[str, str | Path],
    outputs: dict[str, Path | str],
    job_attrs: dict[str, str],
) -> list[BashJob]:
    """
    Run somalier relate on all fingerprints for a single participant,
    then check self-relatedness, alert via Slack, and register results in metamist.

    Returns [relate_job, check_job].
    """
    batch_instance = hail_batch.get_batch()

    # Job 1: somalier relate
    relate_j = batch_instance.new_bash_job(
        f'Somalier identity check {participant_id}',
        job_attrs | {'tool': 'somalier'},
    )
    relate_j.image(config.config_retrieve(['images', 'somalier']))
    storage_gb = 1 + len(somalier_paths) // 4000
    relate_j.storage(f'{storage_gb}GiB')

    relate_j.command('mkdir -p inputs/')
    for sg_id, somalier_path in somalier_paths.items():
        somalier_file = batch_instance.read_input(somalier_path)
        relate_j.command(f'mv {somalier_file} inputs/{sg_id}.somalier')

    relate_j.declare_resource_group(
        output={
            'pairs.tsv': '{root}.pairs.tsv',
            'samples.tsv': '{root}.samples.tsv',
        }
    )
    relate_j.command(f'somalier relate -o {relate_j.output} inputs/*.somalier')
    relate_j.command(f'mv {relate_j.output}.html {relate_j.html_out}')
    batch_instance.write_output(relate_j.output, outputs[f'{participant_id}_prefix'])
    batch_instance.write_output(relate_j.html_out, outputs[f'{participant_id}_html'])

    sg_ids_str = ','.join(sorted(somalier_paths.keys()))

    check_j = batch_instance.new_bash_job(
        f'Somalier identity alert {participant_id}',
        job_attrs,
    )
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))
    check_j.depends_on(relate_j)

    out_html_url = _get_out_html_url(dataset_name, outputs[f'{participant_id}_html'])

    check_j.command(f"""\
python3 -m rd_qc.scripts.check_self_relatedness \\
    --dataset {dataset_name} \\
    --participant-id {participant_id} \\
    --sg-ids {sg_ids_str} \\
    --somalier-pairs {relate_j.output['pairs.tsv']} \\
    --somalier-samples {relate_j.output['samples.tsv']} \\
    --output-pairs {outputs[f'{participant_id}_pairs_tsv']!s} \\
    --output-samples {outputs[f'{participant_id}_samples_tsv']!s} \\
    --output-html {outputs[f'{participant_id}_html']!s} \\
    --html-url {out_html_url} \\
    --output-json {outputs[f'{participant_id}_json']!s}
""")

    return [relate_j, check_j]


def pedigree_check_jobs(  # noqa: PLR0917
    somalier_paths: dict[str, str | Path],
    somalier_self_relatedness_json_paths: list[Path | str],
    outputs: dict[str, Path],
    tmp_prefix: Path,
    dataset_name: str,
    label: str,
    job_attrs: dict[str, str],
) -> list[BashJob]:
    """
    Run somalier relate across all SGs in the dataset with a PED file,
    then validate pedigree, alert via Slack, and register results in metamist.

    Returns [relate_job, check_job].
    """
    batch_instance = hail_batch.get_batch()

    ped_path = outputs['expected_ped']

    # Job 1: somalier relate with --ped --infer
    relate_j = batch_instance.new_bash_job(
        f'Somalier pedigree relate {label}',
        job_attrs | {'tool': 'somalier'},
    )
    relate_j.image(config.config_retrieve(['images', 'somalier']))
    storage_gb = 1 + len(somalier_paths) // 4000
    relate_j.storage(f'{storage_gb}GiB')

    relate_j.command('mkdir -p inputs/')
    for sg_id, somalier_path in somalier_paths.items():
        somalier_file = batch_instance.read_input(somalier_path)
        relate_j.command(f'mv {somalier_file} inputs/{sg_id}.somalier')

    ped_input = batch_instance.read_input(ped_path)

    relate_j.declare_resource_group(
        output={
            'pairs.tsv': '{root}.pairs.tsv',
            'samples.tsv': '{root}.samples.tsv',
        }
    )
    relate_j.command(f'somalier relate --ped {ped_input} -o {relate_j.output} --infer inputs/*.somalier')
    relate_j.command(f'mv {relate_j.output}.html {relate_j.html_out}')
    batch_instance.write_output(relate_j.output, outputs['samples'].parent)
    # First copy of the HTML report written to a uniquely namespaced URL based on this run's AR GUID
    batch_instance.write_output(relate_j.html_out, outputs['html'])
    # Second copy of the HTML report written to the fixed URL
    batch_instance.write_output(relate_j.html_out, outputs['base_html_url'])

    sg_ids_str = ' '.join(sorted(somalier_paths.keys()))
    title = f'Pedigree check [{label}]'

    check_j = batch_instance.new_bash_job(title, job_attrs)
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))
    check_j.depends_on(relate_j)

    out_html_url = _get_out_html_url(dataset_name, outputs['html'])

    cmd = f"""\
python3 -m rd_qc.scripts.check_pedigree \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --sg-ids {sg_ids_str} \\
    --expected-ped {ped_input} \\
    --somalier-samples {relate_j.output['samples.tsv']} \\
    --somalier-pairs {relate_j.output['pairs.tsv']} \\
    --output-pairs {outputs['pairs']!s} \\
    --output-samples {outputs['samples']!s} \\
    --output-html {outputs['html']!s} \\
    --base-output-html {outputs['base_html_url']!s} \\
    --html-url {out_html_url} \\
    --output-json {outputs['json']!s}
touch {check_j.output}
"""
    check_j.command(cmd)
    batch_instance.write_output(check_j.output, outputs['checks'])

    record_j = record_somalier_flags_job(
        dataset_name=dataset_name,
        tmp_prefix=tmp_prefix,
        sg_ids=sg_ids_str,
        somalier_self_relatedness_json_paths=somalier_self_relatedness_json_paths,
        somalier_relatedness_json=str(outputs['json']),
        job_attrs=job_attrs,
    )
    record_j.depends_on(check_j)

    return [relate_j, check_j, record_j]


def record_somalier_flags_job(  # noqa: PLR0917
    dataset_name: str,
    sg_ids: str,
    tmp_prefix: Path,
    somalier_self_relatedness_json_paths: list[Path | str],
    somalier_relatedness_json: str,
    job_attrs: dict | None = None,
) -> BashJob:
    """
    Run job that records all Somalier flags in Metamist by reading the self-relatedness JSON files for
    each sequencing group and the relatedness JSON file for the dataset

    Updates SG meta with any new or changed flags.
    """
    batch_instance = hail_batch.get_batch()
    record_j = batch_instance.new_job('Record Somalier flags', (job_attrs or {}) | {'tool': 'python'})

    record_j.image(config.config_retrieve(['workflow', 'driver_image']))

    # Read in all the self-relatedness JSON files for the SGs
    file_list_path = tmp_prefix / f'{dataset_name}_somalier-self-relatedness-file-list.txt'
    with file_list_path.open('w') as f:
        f.writelines([f'{p}\n' for p in somalier_self_relatedness_json_paths])
    somalier_self_relatedness_jsons = batch_instance.read_input(file_list_path)

    # Read in the full relatedness JSON file for the dataset
    somalier_relatedness_json = batch_instance.read_input(somalier_relatedness_json)

    cmd = f"""\
    mkdir -p self_relatedness_jsons
    cat {somalier_self_relatedness_jsons} | gcloud storage cp -I self_relatedness_jsons/

    python3 -m rd_qc.scripts.record_somalier_flags \\
    --dataset {dataset_name} \\
    --sg-ids {sg_ids} \\
    --somalier-self-relatedness-json-dir self_relatedness_jsons \\
    --somalier-relatedness-json-path {somalier_relatedness_json}
    """

    record_j.command(cmd)
    return record_j
