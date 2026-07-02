"""
Jobs for somalier relate — used by both identity checks and pedigree checks.
"""

from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob


def identity_check_jobs(
    participant_id: str,
    somalier_paths: dict[str, str | Path],
    output_prefix: Path,
    dataset_name: str,
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
    relate_j.storage(f'{storage_gb}Gi')

    cmd = 'mkdir -p inputs/\n'
    for sg_id, somalier_path in somalier_paths.items():
        somalier_file = batch_instance.read_input(str(somalier_path))
        cmd += f'cp {somalier_file} inputs/{sg_id}.somalier\n'

    cmd += f"""\
somalier relate -o results inputs/*.somalier
mv results.pairs.tsv {relate_j.pairs_tsv}
mv results.samples.tsv {relate_j.samples_tsv}
mv results.html {relate_j.html}
"""
    relate_j.command(cmd)

    pairs_out = str(output_prefix) + '.pairs.tsv'
    samples_out = str(output_prefix) + '.samples.tsv'
    html_out = str(output_prefix) + '.html'

    batch_instance.write_output(relate_j.pairs_tsv, pairs_out)
    batch_instance.write_output(relate_j.samples_tsv, samples_out)
    batch_instance.write_output(relate_j.html, html_out)

    # Job 2: check self-relatedness + Slack alert + register in metamist
    kinship_threshold = config.config_retrieve(
        ['workflow', 'somalier_self_check', 'kinship_threshold'],
        0.9,
    )
    sg_ids_str = ','.join(sorted(somalier_paths.keys()))

    check_j = batch_instance.new_bash_job(
        f'Somalier identity alert {participant_id}',
        job_attrs,
    )
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))
    check_j.depends_on(relate_j)

    hail_batch.copy_common_env(check_j)
    hail_batch.authenticate_cloud_credentials_in_job(check_j)
    check_j.command(f"""\
python3 -m rd_qc.scripts.check_self_relatedness \\
    --pairs-tsv {relate_j.pairs_tsv} \\
    --participant-id {participant_id} \\
    --dataset {dataset_name} \\
    --kinship-threshold {kinship_threshold} \\
    --sg-ids {sg_ids_str} \\
    --output-pairs {pairs_out} \\
    --output-samples {samples_out} \\
    --output-html {html_out}
""")

    return [relate_j, check_j]


def pedigree_check_jobs(
    somalier_paths: dict[str, str | Path],
    outputs: dict[str, Path],
    out_html_url: str,
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
    relate_j.storage(f'{storage_gb}Gi')

    cmd = 'mkdir -p inputs/\n'
    for sg_id, somalier_path in somalier_paths.items():
        somalier_file = batch_instance.read_input(str(somalier_path))
        cmd += f'cp {somalier_file} inputs/{sg_id}.somalier\n'

    ped_input = batch_instance.read_input(str(ped_path))
    cmd += f"""\
somalier relate \\
    --ped {ped_input} \\
    -o results \\
    --infer \\
    inputs/*.somalier
mv results.pairs.tsv {relate_j.output_pairs}
mv results.samples.tsv {relate_j.output_samples}
mv results.html {relate_j.output_html}
"""
    relate_j.command(cmd)

    batch_instance.write_output(relate_j.output_pairs, str(outputs['pairs']))
    batch_instance.write_output(relate_j.output_samples, str(outputs['samples']))
    batch_instance.write_output(relate_j.output_html, str(outputs['html']))

    # Job 2: check pedigree + Slack alert + register in metamist
    sg_ids_str = ','.join(sorted(somalier_paths.keys()))
    title = f'Pedigree check [{label}]'

    check_j = batch_instance.new_bash_job(title, job_attrs)
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))
    check_j.depends_on(relate_j)

    hail_batch.copy_common_env(check_j)
    hail_batch.authenticate_cloud_credentials_in_job(check_j)

    cmd = f"""\
python3 -m rd_qc.scripts.check_pedigree \\
    --somalier-samples {relate_j.output_samples} \\
    --somalier-pairs {relate_j.output_pairs} \\
    --ped {ped_input} \\
    --html-url {out_html_url} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --sg-ids {sg_ids_str} \\
    --output-pairs {str(outputs['pairs'])} \\
    --output-samples {str(outputs['samples'])} \\
    --output-html {str(outputs['html'])}
touch {check_j.output}
"""
    check_j.command(cmd)
    batch_instance.write_output(check_j.output, str(outputs['checks']))

    return [relate_j, check_j]