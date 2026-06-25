"""
A job should contain the logic for a single Stage
"""

from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def echo_statement_to_file(statement: str, output_file: Path) -> 'Job':
    """
    This is a simple example of a job that writes a statement to a file.

    Args:
        statement (str): the intended file contents
        output_file (Path): the path to write the file to

    Returns:
        the resulting job
    """
    batch_instance = hail_batch.get_batch()

    # create a job
    j = batch_instance.new_job(f'echo "{statement}" to {output_file}')

    # choose an image to run this job in (default is bare ubuntu)
    j.image(config.config_retrieve(['workflow', 'driver_image']))

    # write the statement to the file
    j.command(f'echo "{statement}" > {j.output}')

    # write the output to the expected location
    batch_instance.write_output(j.output, output_file)

    # return the job
    return j
