"""
Batch jobs to estimate long-read CRAM contamination with VerifyBamID2, from a samtools pileup.
"""

from google.api_core.exceptions import NotFound
from hailtop.batch.job import Job

from rd_qc.utils import get_gcs_object_size

from cpg_flow.filetypes import CramPath
from cpg_utils import Path, config, hail_batch

VERIFYBAMID_THREADS = 4
BASE_STORAGE = 20


def verifybamid_jobs(sg_id: str, cram: CramPath, output: Path, job_attrs: dict) -> list[Job]:
    """
    Pile up the CRAM at the SVD resource's marker sites with samtools, then run VerifyBamID2 on that pileup.

    VerifyBamID2 reading the CRAM itself (--BamFile) was OOM-killed on long reads within the first few Mb of chr1.
    With --PileupFile, samtools reads the long reads and VerifyBamID only parses one short line per marker site.
    Same 1000G SVD resource and --NumPC 4 as align-genotype's short-read CRAM QC, so FREEMIX is comparable.

    Two jobs because the verifybamid image has no samtools. Writes a *.selfSM file, which MultiQC's verifybamid
    module reads for FREEMIX.
    """
    batch_instance = hail_batch.get_batch()
    ref = hail_batch.fasta_res_group(batch_instance)
    svd_prefix = config.config_retrieve(['references', 'contam_svd_prefix'])

    pileup_j = batch_instance.new_bash_job(f'Long-read VerifyBamID pileup {sg_id}', job_attrs | {'tool': 'samtools'})
    pileup_j.image(config.config_retrieve(['images', 'samtools']))
    try:
        storage_gb = get_gcs_object_size(str(cram.path)) + BASE_STORAGE
    except NotFound:
        storage_gb = 100
    pileup_j.storage(f'{storage_gb}GiB')

    marker_bed = batch_instance.read_input(f'{svd_prefix}.bed')
    cram_localised = batch_instance.read_input_group(
        cram=str(cram.path),
        crai=str(cram.index_path or f'{cram.path}.crai'),
    ).cram

    # -B: no BAQ, which is tuned for Illumina errors and penalises bases around long-read indel errors
    # -q 20: drop ambiguously mapped reads, so paralogous regions don't look like contamination
    # --positions: only the SVD resource's own marker sites, so the pileup lines up with what VerifyBamID models
    pileup_j.command(f"""\
    samtools mpileup -B -q 20 -f {ref.base} --positions {marker_bed} -o {pileup_j.pileup} {cram_localised}
    """)

    verify_j = batch_instance.new_bash_job(f'Long-read VerifyBamID {sg_id}', job_attrs | {'tool': 'VerifyBamID'})
    verify_j.image(config.config_retrieve(['images', 'verifybamid']))
    verify_j.cpu(VERIFYBAMID_THREADS)

    verify_ref = hail_batch.fasta_res_group(batch_instance)
    svd = batch_instance.read_input_group(UD=f'{svd_prefix}.UD', mu=f'{svd_prefix}.mu', bed=f'{svd_prefix}.bed')

    # The SVD files are localised side by side under their original names, so the prefix is the UD path minus .UD
    verify_j.command(f"""\
    SVD_PREFIX=$(echo {svd.UD} | sed 's/\\.UD$//')
    VerifyBamID \\
        --NumThread {VERIFYBAMID_THREADS} \\
        --Verbose \\
        --NumPC 4 \\
        --Output OUTPUT \\
        --PileupFile {pileup_j.pileup} \\
        --SVDPrefix $SVD_PREFIX \\
        --Reference {verify_ref.base} \\
        1>/dev/null

    cp OUTPUT.selfSM {verify_j.out_selfsm}
    """)

    batch_instance.write_output(verify_j.out_selfsm, str(output))
    return [pileup_j, verify_j]
