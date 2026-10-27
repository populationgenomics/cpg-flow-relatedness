"""
Batch job to estimate long-read CRAM contamination with VerifyBamID2.
"""

from google.api_core.exceptions import NotFound
from hailtop.batch.job import Job

from rd_qc.utils import get_gcs_object_size

from cpg_flow.filetypes import CramPath
from cpg_utils import Path, config, hail_batch

VERIFYBAMID_THREADS = 4


def verifybamid_job(sg_id: str, cram: CramPath, output: Path, job_attrs: dict) -> Job:
    """
    Run VerifyBamID2 directly on the CRAM, as align-genotype's short-read CRAM QC does, with the same 1000G SVD
    resource and --NumPC 4, so long-read FREEMIX is measured the same way.

    Writes a *.selfSM file, which MultiQC's verifybamid module reads for FREEMIX.
    """
    batch_instance = hail_batch.get_batch()

    j = batch_instance.new_bash_job(f'Long-read VerifyBamID {sg_id}', job_attrs | {'tool': 'VerifyBamID'})
    j.image(config.config_retrieve(['images', 'verifybamid']))
    j.cpu(VERIFYBAMID_THREADS)
    try:
        storage_gb = get_gcs_object_size(str(cram.path))
    except NotFound:
        storage_gb = 100
    j.storage(f'{storage_gb}GiB')

    ref = hail_batch.fasta_res_group(batch_instance)
    svd_prefix = config.config_retrieve(['references', 'contam_svd_prefix'])
    svd = batch_instance.read_input_group(ud=f'{svd_prefix}.UD', mu=f'{svd_prefix}.mu', bed=f'{svd_prefix}.bed')
    cram_localised = batch_instance.read_input_group(
        cram=str(cram.path),
        crai=str(cram.index_path or f'{cram.path}.crai'),
    ).cram

    j.command(f"""\
    VerifyBamID \\
        --NumThread {VERIFYBAMID_THREADS} \\
        --Verbose \\
        --NumPC 4 \\
        --Output OUTPUT \\
        --BamFile {cram_localised} \\
        --Reference {ref.base} \\
        --UDPath {svd.ud} \\
        --MeanPath {svd.mu} \\
        --BedPath {svd.bed} \\
        1>/dev/null

    cp OUTPUT.selfSM {j.out_selfsm}
    """)

    batch_instance.write_output(j.out_selfsm, str(output))
    return j
