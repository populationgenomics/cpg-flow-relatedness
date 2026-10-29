"""
Stages for the long-read rd_qc workflow: long-read CRAM QC, and the somalier pedigree check and flags report
re-declared to run after it.

Run with run_workflow_lrs, alongside a config override setting workflow.sequencing_technology = 'long-read' and the
long-read input_cohorts. Short-read CRAM QC happens upstream, in align-genotype, so the short-read workflow
(run_workflow) has no equivalent of these.
"""

from rd_qc.jobs import lrs_multiqc, lrs_samtools_stats, lrs_verifybamid
from rd_qc.stages import (
    GenerateMissingSomalierFingerprints,
    GenerateSomalierFlagsReport,
    SomalierPedigreeCheck,
    SomalierSelfCheck,
)
from rd_qc.utils import convert_to_web_url

from cpg_flow import stage, targets
from cpg_utils import Path, to_path


@stage.stage(analysis_type='qc', analysis_keys=['stats'])
class LongReadQcSamtoolsStats(stage.SequencingGroupStage):
    """
    samtools stats on a representative sample of a long-read CRAM, for error_rate, plus whole-CRAM idxstats for
    reads_mapped_percent.

    The registered `qc` analysis carries meta.stage, which LongReadMultiQC queries on.
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        prefix = sequencing_group.dataset.prefix() / 'qc' / 'long_read' / 'samtools_stats' / sequencing_group.id
        return {
            'stats': to_path(f'{prefix}.samtools-stats'),
            'idxstats': to_path(f'{prefix}.idxstats'),
            'mapped_mqc': to_path(f'{prefix}.lrs_mapped_mqc.tsv'),
        }

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        outputs = self.expected_outputs(sequencing_group)
        if sequencing_group.cram is None:
            return self.make_outputs(sequencing_group, error_msg=f'{sequencing_group.id}: no CRAM in Metamist')

        job = lrs_samtools_stats.samtools_stats_job(
            sg_id=sequencing_group.id,
            cram=sequencing_group.cram,
            outputs=outputs,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        return self.make_outputs(sequencing_group, data=outputs, jobs=[job], meta={'stage': self.name})


@stage.stage(analysis_type='qc')
class LongReadQcVerifyBamId(stage.SequencingGroupStage):
    """
    VerifyBamID2 on a long-read CRAM, for FREEMIX.

    The registered `qc` analysis carries meta.stage, which LongReadMultiQC queries on.
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> Path:
        return (
            sequencing_group.dataset.prefix()
            / 'qc'
            / 'long_read'
            / 'verifybamid'
            / f'{sequencing_group.id}.verify-bamid.selfSM'
        )

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        output = self.expected_outputs(sequencing_group)
        if sequencing_group.cram is None:
            return self.make_outputs(sequencing_group, error_msg=f'{sequencing_group.id}: no CRAM in Metamist')

        job = lrs_verifybamid.verifybamid_job(
            sg_id=sequencing_group.id,
            cram=sequencing_group.cram,
            output=output,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        return self.make_outputs(sequencing_group, data=output, jobs=[job], meta={'stage': self.name})


@stage.stage(
    required_stages=[LongReadQcSamtoolsStats, LongReadQcVerifyBamId],
    analysis_type='qc',
    analysis_keys=['json', 'html'],
)
class LongReadMultiQC(stage.DatasetStage):
    """
    MultiQC over the dataset's long-read CRAM QC: this run's samtools stats and VerifyBamID outputs, plus earlier
    runs' and the ingested mosdepth outputs from Metamist. Checks them against the per-platform thresholds and
    records flags in SG meta under `cram_qc_flags`.
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        sg_hash = dataset.get_alignment_inputs_hash()
        return {
            'latest': dataset.web_prefix() / 'qc' / 'long_read' / 'latest' / 'multiqc.html',
            'html': dataset.web_prefix() / 'qc' / 'long_read' / sg_hash / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'long_read' / sg_hash / 'multiqc_data.json',
            'checks': dataset.prefix() / 'qc' / 'long_read' / sg_hash / 'qc-checks.json',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        # This run's outputs, which Metamist won't have registered yet when lrs_multiqc queries it
        current_paths: dict[str, list[str]] = {}
        for sg in dataset.get_sequencing_groups():
            stats = inputs.as_dict(sg, LongReadQcSamtoolsStats)
            current_paths[sg.id] = [
                str(stats['stats']),
                str(stats['mapped_mqc']),
                str(inputs.as_path(sg, LongReadQcVerifyBamId)),
            ]

        jobs = lrs_multiqc.lrs_multiqc(
            dataset=dataset,
            outputs=outputs,
            out_html_url=convert_to_web_url(dataset.name, outputs['html']),
            label='CRAM',
            job_attrs=self.get_job_attrs(dataset),
            current_paths=current_paths,
        )
        return self.make_outputs(dataset, data=outputs, jobs=jobs)


# The pedigree check and flags report, re-declared for long-read runs. The logic is the short-read stages', reused
# through `__wrapped__`: cpg-flow's @stage decorator returns a factory function, and keeps the class there.


@stage.stage(required_stages=[GenerateMissingSomalierFingerprints, SomalierSelfCheck, LongReadMultiQC])
class LongReadSomalierPedigreeCheck(SomalierPedigreeCheck.__wrapped__):  # type: ignore[name-defined,misc]
    """
    SomalierPedigreeCheck, run after LongReadMultiQC.

    check_pedigree.py reads unresolved FREEMIX flags from Metamist when it runs and excludes contaminated SGs from
    relatedness flags. Long-read FREEMIX comes from LongReadMultiQC in the same run, so this waits for it.
    """


@stage.stage(required_stages=[LongReadSomalierPedigreeCheck], forced=True)
class LongReadSomalierFlagsReport(GenerateSomalierFlagsReport.__wrapped__):  # type: ignore[name-defined,misc]
    """GenerateSomalierFlagsReport, run after the long-read pedigree check."""
