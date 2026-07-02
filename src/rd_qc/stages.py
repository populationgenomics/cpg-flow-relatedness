"""Stages for the rd_qc somalier QC workflow."""

from argparse import ArgumentParser

from cpg_flow import stage, targets, workflow
from cpg_utils import config, to_path, Path

from rd_qc.jobs import generate_somalier, relate
from rd_qc.utils import (
    build_ped_content,
    find_sgids_without_somalier,
    get_project_sgs_and_fingerprints,
    select_somalier_extract_targets,
    sg_ids_tag,
)


@stage.stage()
class GenerateMissingSomalierFingerprints(stage.DatasetStage):

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        index = get_project_sgs_and_fingerprints(dataset.name)
        missing_sgids = find_sgids_without_somalier(index)

        if not missing_sgids:
            return {}

        extract_targets = select_somalier_extract_targets(
            dataset.name,
            tuple(sorted(missing_sgids)),
        )

        return {sg_id: to_path(f'{source_file}.somalier') for sg_id, source_file in extract_targets.items()}

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        if not outputs:
            return self.make_outputs(dataset, data=outputs)

        index = get_project_sgs_and_fingerprints(dataset.name)
        missing_sgids = find_sgids_without_somalier(index)
        extract_targets = select_somalier_extract_targets(
            dataset.name,
            tuple(sorted(missing_sgids)),
        )

        jobs = generate_somalier.somalier_jobs(
            somalier_targets=extract_targets,
            somalier_outputs=outputs,
            project=dataset.name,
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage.stage(required_stages=[GenerateMissingSomalierFingerprints])
class RunCrossTypeIdentityChecks(stage.DatasetStage):

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        index = get_project_sgs_and_fingerprints(dataset.name)
        output_prefix = to_path(config.config_retrieve(['storage', dataset.name, 'default'])) / 'identity_checks'

        outputs = {}
        for participant_id, sg_list in index.by_participant.items():
            if len(sg_list) < 2:
                continue

            tag = sg_ids_tag([info.sg_id for info in sg_list])
            prefix = output_prefix / participant_id / f'{tag}.somalier_identity_check'
            outputs[f'{participant_id}_pairs_tsv'] = to_path(str(prefix) + '.pairs.tsv')
            outputs[f'{participant_id}_samples_tsv'] = to_path(str(prefix) + '.samples.tsv')
            outputs[f'{participant_id}_html'] = to_path(str(prefix) + '.html')

        return outputs

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        if not outputs:
            return self.make_outputs(dataset, data=outputs)

        index = get_project_sgs_and_fingerprints(dataset.name)

        # Patch in newly-generated fingerprints via O(1) sg_id lookup
        new_fingerprints = inputs.as_dict(dataset, GenerateMissingSomalierFingerprints)
        for sg_id, new_path in new_fingerprints.items():
            if sg_id in index.by_sg:
                index.by_sg[sg_id].somalier_path = str(new_path)

        output_prefix = to_path(config.config_retrieve(['storage', dataset.name, 'default'])) / 'identity_checks'

        all_jobs = []
        for participant_id, sg_list in index.by_participant.items():
            if len(sg_list) < 2:
                continue

            somalier_paths = {info.sg_id: info.somalier_path for info in sg_list if info.somalier_path is not None}
            if len(somalier_paths) < 2:
                continue

            tag = sg_ids_tag([info.sg_id for info in sg_list])
            prefix = output_prefix / participant_id / f'{tag}.somalier_identity_check'

            jobs = relate.identity_check_jobs(
                participant_id=participant_id,
                somalier_paths=somalier_paths,
                output_prefix=prefix,
                dataset_name=dataset.name,
                job_attrs={'participant': participant_id},
            )
            all_jobs.extend(jobs)

        return self.make_outputs(dataset, data=outputs, jobs=all_jobs)


@stage.stage(required_stages=[GenerateMissingSomalierFingerprints])
class SomalierPedigreeCheck(stage.DatasetStage):

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        prefix = to_path(config.config_retrieve(['storage', dataset.name, 'default'])) / 'somalier_checks' / 'pedigree'
        web_prefix = dataset.web_prefix() / 'somalier_checks' / 'pedigree'

        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            'expected_ped': prefix / f'{dataset.name}.expected.ped',
            'html': web_prefix / 'somalier-pedigree.html',
            'checks': prefix / f'{dataset.name}-checks.done',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        index = get_project_sgs_and_fingerprints(dataset.name)

        # Patch in newly-generated fingerprints via O(1) sg_id lookup
        new_fingerprints = inputs.as_dict(dataset, GenerateMissingSomalierFingerprints)
        for sg_id, new_path in new_fingerprints.items():
            if sg_id in index.by_sg:
                index.by_sg[sg_id].somalier_path = str(new_path)

        # Collect all somalier paths — participant_id comes from the dataclass
        somalier_paths: dict[str, str] = {}
        for info in index.by_sg.values():
            if info.somalier_path is not None:
                somalier_paths[info.sg_id] = info.somalier_path

        # Build PED file content and write to GCS at orchestration time
        ped_content = build_ped_content(dataset.name, index)
        ped_path = outputs['expected_ped']
        with to_path(ped_path).open('w') as f:
            f.write(ped_content)

        html_url = str(outputs['html']).replace(
            str(dataset.web_prefix()),
            dataset.web_url(),
        )

        jobs = relate.pedigree_check_jobs(
            somalier_paths=somalier_paths,
            outputs=outputs,
            out_html_url=html_url,
            dataset_name=dataset.name,
            label=f'{dataset.name} Somalier',
            job_attrs={},
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)