"""This file exists to define all the Stages for the workflow."""

from itertools import combinations

from cpg_flow import stage, targets, workflow, utils as flow_utils
from cpg_flow.targets import dataset
from cpg_utils import config, to_path, Path

from rd_qc.jobs import generate_somalier
from rd_qc.utils import get_project_sgs_and_fingerprints, find_sgids_without_somalier, select_somalier_extract_targets


@stage.stage()
class GenerateMissingSomalierFingerprints(stage.DatasetStage):

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        """
        Ok, right, this is where it gets weird. First stage, weird already...
        """

        # get all SG IDs in the dataset. Not in the 'Dataset', but across all sequencing groups in the whole project
        all_sgid_somaliers = get_project_sgs_and_fingerprints(dataset.name)

        all_sgids_missing_somalier = find_sgids_without_somalier(...)

        new_somalier_targets = select_somalier_extract_targets(...)

        return {key: to_path(f'{value}.somalier') for key, value in new_somalier_targets.items()}

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        """
        This is where we generate jobs for this stage.
        """

        # all the required outputs. Inputs are this - the somalier extension
        outputs = self.expected_outputs(dataset)

        # or just re-generate the inputs using the cached methods...
        all_sgid_somaliers = get_project_sgs_and_fingerprints(dataset.name)

        all_sgids_missing_somalier = find_sgids_without_somalier(...)

        new_somalier_targets = select_somalier_extract_targets(...)

        # feed the new_somalier_targets into a job to generate each file

        # run a second job to register each somalier file (just like in the rna dashboard)
        jobs = generate_somalier.somalier_jobs(new_somalier_targets)


        # return the jobs and outputs
        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage.stage(required_stages=[GenerateMissingSomalierFingerprints])
class RunCrossTypeIdentityChecks(stage.DatasetStage):

    def expected_outputs(self, dataset: targets.Dataset) -> list[Path]:
        """
        For this dataset we first need to know which pairs are possible. We do this by re-using the cached db query
        This contains all possible combinations, and the expectation that somalier files that didn't exist were populated by the previous stage.
        """
        all_sgid_somaliers = get_project_sgs_and_fingerprints(dataset.name)
        output_prefix = to_path(config.config_retrieve(['storage', dataset.name, 'default'])) / 'identity_checks'

        output_paths = []
        for participant_id, sgid_map in all_sgid_somaliers.items():

            for sgid1, sgid2 in combinations(sgid_map.keys(), 2):

                # re-use the sgid sorting method to get a consistent file name
                new_file = output_prefix / participant_id / f'{sgid1}_{sgid2}.check'

                # check if this already exists using the efficient cached existence method
                if not flow_utils.exists(new_file):
                    # if it doesn't exist, add it as an output
                    output_paths.append(new_file)

        return output_paths

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        """
        This might be a bit rogue
        """

        outputs = self.expected_outputs(dataset)

        # get just the file names as a nifty little set
        output_names = {out.name for out in outputs}

        # take the input from the previous stage ({sgid: new somalier file})
        new_somalier_files = inputs.as_dict(dataset, GenerateMissingSomalierFingerprints)

        all_sgid_somaliers = get_project_sgs_and_fingerprints(dataset.name)

        # update the None entries in all_sgid_somaliers with the new somalier fingerprints generated in the previous stage
        ...

        for participant, sgid_map in all_sgid_somaliers.items():
            for sgid1, sgid2 in combinations(sgid_map.keys(), 2):
                # if this one doesn't exist yet (in list of non-existent files)
                if f'{sgid1}_{sgid2}.check' in output_names:

                    # invoke the method which takes two somalier files, runs _relate_, creates output, and registers
                    ...

