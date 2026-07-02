
"""
Main entry point for the rd_qc workflow.
Imports all stages and begins the CPG-Flow stage discovery and graph construction process.
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow

from rd_qc.stages import (
    GenerateMissingSomalierFingerprints,
    RunCrossTypeIdentityChecks,
    SomalierPedigreeCheck,
)


def cli_main() -> None:
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    stages = [GenerateMissingSomalierFingerprints, RunCrossTypeIdentityChecks, SomalierPedigreeCheck]
    run_workflow(name='rd_qc', stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()