"""
Batch job to collect long-read alignment stats with samtools: stats on a representative sample, plus
whole-CRAM mapped/unmapped counts.
"""

from google.api_core.exceptions import NotFound
from hailtop.batch.job import Job

from rd_qc.utils import get_gcs_object_size

from cpg_flow.filetypes import CramPath
from cpg_utils import Path, config, hail_batch

SAMTOOLS_THREADS = 4

# Evenly spaced windows across chr1-22: the autosomes are laid end to end, n windows are placed at equal steps along
# that total, and each is clipped to its chromosome's end.
WINDOWS_AWK = r"""awk -v n={n_windows} -v w={window_bp} 'BEGIN {{ OFS = "\t" }}
$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {{ name[++c] = $1; len[c] = $2; total += $2 }}
END {{
    step = total / n; pos = step / 2; offset = 0
    for (i = 1; i <= c; i++) {{
        while (pos < offset + len[i]) {{
            start = int(pos - offset); end = start + w
            if (end > len[i]) end = len[i]
            print name[i], start, end
            pos += step
        }}
        offset += len[i]
    }}
}}' {fai} > windows.bed"""

# reads_mapped_percent from whole-CRAM idxstats, written as MultiQC custom content so it reaches the general stats
# table. The windowed sample only ever sees mapped reads, so its own reads_mapped_percent is always ~100%.
MAPPED_MQC_HEADER = """\
# id: 'lrs_idxstats'
# section_name: 'Long-read mapping (whole CRAM, samtools idxstats)'
# plot_type: 'generalstats'
Sample\treads_mapped_percent"""


def samtools_stats_job(sg_id: str, cram: CramPath, outputs: dict[str, Path], job_attrs: dict) -> Job:
    """
    samtools stats on ~1M reads from evenly spaced autosomal windows, read via the CRAM index rather than a full pass,
    for error_rate and the read length distribution. samtools idxstats on the whole CRAM for reads_mapped_percent.
    """
    batch_instance = hail_batch.get_batch()

    j = batch_instance.new_bash_job(f'Long-read samtools stats {sg_id}', job_attrs | {'tool': 'samtools'})
    j.image(config.config_retrieve(['images', 'samtools']))
    j.cpu(SAMTOOLS_THREADS)
    try:
        storage_gb = get_gcs_object_size(str(cram.path))
    except NotFound:
        storage_gb = 100
    j.storage(f'{storage_gb}GiB')

    ref = hail_batch.fasta_res_group(batch_instance)
    cram_localised = batch_instance.read_input_group(
        cram=str(cram.path),
        crai=str(cram.index_path or f'{cram.path}.crai'),
    ).cram

    windows_cmd = WINDOWS_AWK.format(
        n_windows=config.config_retrieve(['lrs_qc', 'samtools_stats', 'n_windows'], 500),
        window_bp=config.config_retrieve(['lrs_qc', 'samtools_stats', 'window_bp'], 1_000_000),
        fai=ref.fai,
    )

    j.command(f"""\
    {windows_cmd}

    samtools view -u -M -L windows.bed --reference {ref.base} {cram_localised} \\
        | samtools stats -@ {SAMTOOLS_THREADS - 1} - > {j.stats}

    samtools idxstats {cram_localised} > {j.idxstats}

    cat > {j.mapped_mqc} <<'EOF'
{MAPPED_MQC_HEADER}
EOF
    awk -v sg={sg_id} '{{ mapped += $3; unmapped += $4 }}
        END {{ printf "%s\\t%.4f\\n", sg, (mapped + unmapped) ? 100 * mapped / (mapped + unmapped) : 0 }}' \\
        {j.idxstats} >> {j.mapped_mqc}
    """)

    batch_instance.write_output(j.stats, str(outputs['stats']))
    batch_instance.write_output(j.idxstats, str(outputs['idxstats']))
    batch_instance.write_output(j.mapped_mqc, str(outputs['mapped_mqc']))
    return j
