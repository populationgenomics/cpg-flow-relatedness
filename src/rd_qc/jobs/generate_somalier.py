from hailtop.batch.job import BashJob

from cpg_flow.status import complete_analysis_job
from cpg_utils import to_path, Path

"""
Job to generate somalier fingerprints for SGs missing them.
Runs somalier extract on each input file and registers the result in metamist.
"""

import re

from cpg_flow.status import complete_analysis_job
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob

_CRAM_PATTERN = re.compile(r'\.cram$')


def somalier_jobs(
    somalier_targets: dict[str, str],
    somalier_outputs: dict[str, Path],
    project: str,
) -> list[BashJob]:
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

        j = batch_instance.new_bash_job(
            f'Somalier extract {sg_id}',
            {'tool': 'somalier', 'sg': sg_id},
        )
        j.image(config.config_retrieve(['images', 'somalier']))

        #crams need more space
        is_cram = bool(_CRAM_PATTERN.search(source_file))
        if is_cram:
            storage_gb = config.config_retrieve(
                ['workflow', 'resource_overrides', 'somalier_extract', 'storage_gib'],
                50,
            )
            j.storage(f'{storage_gb}GB')
            localised = batch_instance.read_input_group(
                cram=source_file,
                crai=f'{source_file}.crai',
            ).cram
        else:
            j.storage('10GB')
            localised = batch_instance.read_input_group(
                **{'vcf.gz': source_file, 'vcf.gz.tbi': f'{source_file}.tbi'},
            )['vcf.gz']

        j.command(f"""\
        export SOMALIER_SAMPLE_NAME={sg_id}
        somalier extract -d extracted/ --sites {sites} -f {ref.base} {localised}
        mv extracted/*.somalier {j.output_file}
        """)

        batch_instance.write_output(j.output_file, str(output_path))

        complete_analysis_job(
            batch=batch_instance,
            output=str(output_path),
            sequencing_group_ids=[sg_id],
            analysis_type='somalier',
            meta={},
            depends_on=j,
            project_name=project,
        )

        jobs.append(j)

    return jobs

#def somalier_jobs(somalier_targets: dict[str, Path]) -> list[BashJob]:
#    """
#    Method to take all the somalier targets, decide what type of index each has (tbi, crai) based on file type.
#    localises the main file and index inside a new somalier job, runs somalier extract
#    writes the new file to f'{input}.somalier'
#    registers the result using complete_analysis_job, writing a single-SG ID 'somalier' analysis type
#    """
