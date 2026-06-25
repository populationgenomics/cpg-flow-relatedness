from hailtop.batch.job import BashJob

from cpg_flow.status import complete_analysis_job
from cpg_utils import to_path, Path


def somalier_jobs(somalier_targets: dict[str, Path]) -> list[BashJob]:
    """
    Method to take all the somalier targets, decide what type of index each has (tbi, crai) based on file type.
    localises the main file and index inside a new somalier job, runs somalier extract
    writes the new file to f'{input}.somalier'
    registers the result using complete_analysis_job, writing a single-SG ID 'somalier' analysis type
    """
    ...
