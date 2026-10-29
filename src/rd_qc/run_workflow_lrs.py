"""
Entry point for the long-read rd_qc workflow: long-read CRAM QC, then the somalier checks.

Run with a config override setting workflow.sequencing_technology = 'long-read' and the long-read input_cohorts,
which scopes the pedigree check, the flags report and the input SGs to long-read data.
"""

from argparse import ArgumentParser

from rd_qc.lrs_stages import (
    LongReadMultiQC,
    LongReadQcSamtoolsStats,
    LongReadQcVerifyBamId,
    LongReadSomalierFlagsReport,
    LongReadSomalierPedigreeCheck,
)
from rd_qc.stages import GenerateMissingSomalierFingerprints, SomalierSelfCheck

from cpg_flow.workflow import run_workflow


def cli_main() -> None:
    """
    CLI entrypoint - starts up the long-read workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    stages = [
        GenerateMissingSomalierFingerprints,
        SomalierSelfCheck,
        LongReadQcSamtoolsStats,
        LongReadQcVerifyBamId,
        LongReadMultiQC,
        LongReadSomalierPedigreeCheck,
        LongReadSomalierFlagsReport,
    ]
    run_workflow(name='rd_qc', stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
