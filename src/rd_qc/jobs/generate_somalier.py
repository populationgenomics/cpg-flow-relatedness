"""
Job to generate somalier fingerprints for SGs missing them.
Runs somalier extract on each input file and registers the result in metamist.
"""

from google.api_core.exceptions import NotFound
from hailtop.batch.job import Job

from rd_qc.utils import get_gcs_object_size

from cpg_flow.status import complete_analysis_job
from cpg_utils import Path, config, hail_batch


def register_analyses(output, analysis_type, cohort_ids, sg_ids, project_name, meta):
    complete_analysis_job(output, analysis_type, cohort_ids, sg_ids, project_name, meta)


def somalier_jobs(
    somalier_targets: dict[str, str],
    somalier_outputs: dict[str, Path],
    sg_id_map: dict[str, dict],
    project: str,
) -> list[Job]:
    """
    For each SG needing a fingerprint, run somalier extract and register the result.

    Args:
        somalier_targets: {sg_id: source_file_path} — files to extract from
        somalier_outputs: {sg_id: output .somalier path} — where to write outputs
        project: metamist project name for registration
    """
    batch_instance = hail_batch.get_batch()
    ref = hail_batch.fasta_res_group(batch_instance)
    sites = batch_instance.read_input(config.config_retrieve(['references', 'somalier_sites']))

    jobs = []
    for sg_id, source_file in somalier_targets.items():
        output_path = somalier_outputs[sg_id]
        sample_id = sg_id_map[sg_id]['sample_external_id']
        participant_id = sg_id_map[sg_id]['participant_external_id']
        j = batch_instance.new_bash_job(
            f'{project} Somalier extract {sg_id} | {sample_id} | {participant_id}',
            {'tool': 'somalier', 'sg': sg_id},
        )
        j.image(config.config_retrieve(['images', 'somalier']))
        try:
            storage_gb = get_gcs_object_size(source_file)
        except NotFound:
            storage_gb = 50
        j.storage(f'{storage_gb}GiB')

        is_cram = source_file.endswith('.cram')
        if is_cram:
            localised = batch_instance.read_input_group(
                cram=source_file,
                crai=f'{source_file}.crai',
            ).cram
        else:
            localised = batch_instance.read_input_group(
                **{'vcf.gz': source_file, 'vcf.gz.tbi': f'{source_file}.tbi'},
            )['vcf.gz']

        j.command(f"""\
        export SOMALIER_SAMPLE_NAME={sg_id}
        somalier extract -d extracted/ --sites {sites} -f {ref.base} {localised}
        mv extracted/*.somalier {j.output_file}
        """)

        batch_instance.write_output(j.output_file, str(output_path))

        registration_job = batch_instance.new_python_job(
            f'Register somalier {sg_id}',
            attributes={'tool': 'metamist'},
        )
        registration_job.image(config.config_retrieve(['workflow', 'driver_image']))

        registration_job.call(
            register_analyses,
            output=str(output_path),
            analysis_type='somalier',
            cohort_ids=[],
            sg_ids=[sg_id],
            project_name=project,
            meta={},
        )
        registration_job.depends_on(j)

        jobs.append(j)
        jobs.append(registration_job)

    return jobs
