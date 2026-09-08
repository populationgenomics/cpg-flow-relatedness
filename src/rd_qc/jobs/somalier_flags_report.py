"""
Batch job to generate the somalier flags HTML report from Metamist somalier flags.
"""

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch


def somalier_flags_report_job(
    dataset: str,
    outputs: dict[str, Path],
    out_html_url: str,
    somalier_html_url: str,
    job_attrs: dict,
) -> Job:
    """
    Create a Hail Batch job that queries Metamist for all relatedness flags in the dataset.
    Generates a summary HTML report showing all open and resolved relatedness flags.
    """
    batch = hail_batch.get_batch()

    j = batch.new_bash_job(f'Somalier Flags Report: {dataset}', job_attrs | {'tool': 'python'})
    j.image(config.config_retrieve(['workflow', 'driver_image'])).memory('standard').cpu(2)

    j.command(
        f"""\
    python3 -m align_genotype.scripts.somalier_flags_report \\
        --dataset {dataset} \\
        --output-html {outputs['timestamped']} \\
        --base-output-html {outputs['html']} \\
        --flags-html-url {out_html_url} \\
        --somalier-html-url {somalier_html_url}
    """
    )

    return j
