"""
trivially simple python job example, using a utility constant
"""

from typing import TYPE_CHECKING

from workflow_name.utils import DATE_STRING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def print_file_contents(input_file: str) -> str:
    """
    This is a simple example of a job that prints the contents of a file.

    Args:
        input_file (str): the path to the file to print
    """

    with open(input_file) as f:
        contents = f.read()

    print(f'Contents of {input_file} on {DATE_STRING}:')
    print(contents)
    return contents


def set_up_printing_python_job(input_file: str, output_file: Path, job_attrs: dict[str, str]) -> 'Job':
    """
    This is a simple example of a job that prints the contents of a file.
    This is the logic for the stage, calling the pythonJob to do the work

    Args:
        input_file (str): the path to the file to print
        output_file (Path): the path to write the result to
        job_attrs (dict[str, str]): attributes to attach to the job
    """

    batch_instance = hail_batch.get_batch()

    # localise the file
    local_input = batch_instance.read_input(input_file)

    # run the PythonJob
    job = batch_instance.new_python_job(f'Read {input_file}', attributes=job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    pyjob_output = job.call(
        print_file_contents,
        local_input,
    )
    batch_instance.write_output(pyjob_output.as_str(), output_file)
    return job
