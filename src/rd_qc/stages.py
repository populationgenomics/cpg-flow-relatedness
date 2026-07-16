"""Stages for the rd_qc somalier QC workflow."""

from rd_qc.jobs import generate_somalier, relate
from rd_qc.utils import (
    SomalierIndex,
    build_ped_content,
    find_sgids_without_somalier,
    get_project_sgs_and_fingerprints,
    select_somalier_extract_targets,
    sg_ids_tag,
)

from cpg_flow import stage, targets
from cpg_utils import Path, to_path
from cpg_utils.existence_checks import exists

_MIN_SGS_FOR_IDENTITY_CHECK = 2


def _relevant_sg_ids(dataset: targets.Dataset, index: SomalierIndex) -> set[str]:
    """SG IDs for participants who have at least one SG in the cohort."""
    cohort_sg_ids = {sg.id for sg in dataset.get_sequencing_groups()}
    cohort_participants = {index.by_sg[sg_id].participant_id for sg_id in cohort_sg_ids if sg_id in index.by_sg}
    return {info.sg_id for pid in cohort_participants for info in index.by_participant[pid]}


@stage.stage()
class GenerateMissingSomalierFingerprints(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        index = get_project_sgs_and_fingerprints(dataset.name)
        outputs: dict[str, Path] = {}

        for info in index.by_sg.values():
            if info.somalier_path is not None:
                outputs[info.sg_id] = to_path(info.somalier_path)

        relevant = _relevant_sg_ids(dataset, index)
        missing_sgids = find_sgids_without_somalier(index) & relevant
        if missing_sgids:
            extract_targets = select_somalier_extract_targets(
                dataset.name,
                tuple(sorted(missing_sgids)),
            )
            for sg_id, source_file in extract_targets.items():
                outputs[sg_id] = to_path(f'{source_file}.somalier')

        return outputs

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        outputs = self.expected_outputs(dataset)

        index = get_project_sgs_and_fingerprints(dataset.name)
        missing_sgids = find_sgids_without_somalier(index) & _relevant_sg_ids(dataset, index)

        if not missing_sgids:
            return self.make_outputs(dataset, data=outputs)

        extract_targets = select_somalier_extract_targets(
            dataset.name,
            tuple(sorted(missing_sgids)),
        )

        # Only pass the missing SGs' output paths to the job builder
        somalier_outputs = {sg_id: outputs[sg_id] for sg_id in extract_targets}

        jobs = generate_somalier.somalier_jobs(
            somalier_targets=extract_targets,
            somalier_outputs=somalier_outputs,
            project=dataset.name,
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage.stage(required_stages=[GenerateMissingSomalierFingerprints])
class RunCrossTypeIdentityChecks(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        index = get_project_sgs_and_fingerprints(dataset.name)
        output_prefix = dataset.prefix() / 'identity_checks'
        web_output_prefix = dataset.web_prefix() / 'identity_checks'

        outputs = {}
        for participant_id, sg_list in index.by_participant.items():
            if len(sg_list) < _MIN_SGS_FOR_IDENTITY_CHECK:
                continue

            tag = sg_ids_tag([info.sg_id for info in sg_list])
            prefix = output_prefix / participant_id / f'{tag}.somalier_identity_check'
            web_prefix = web_output_prefix / participant_id / f'{tag}.somalier_identity_check'
            outputs[f'{participant_id}_prefix'] = str(prefix)
            outputs[f'{participant_id}_pairs_tsv'] = to_path(str(prefix) + '.pairs.tsv')
            outputs[f'{participant_id}_samples_tsv'] = to_path(str(prefix) + '.samples.tsv')
            outputs[f'{participant_id}_html'] = to_path(str(web_prefix) + '.html')

        return outputs

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        if not outputs:
            return self.make_outputs(dataset, data=outputs)

        all_somalier = inputs.as_dict(dataset, GenerateMissingSomalierFingerprints)
        index = get_project_sgs_and_fingerprints(dataset.name)

        missing_participants = {
            pid
            for pid, sg_list in index.by_participant.items()
            if len(sg_list) >= _MIN_SGS_FOR_IDENTITY_CHECK and not exists(outputs[f'{pid}_samples_tsv'])
        }

        if not missing_participants:
            return self.make_outputs(dataset, data=outputs)

        all_jobs = []
        for participant_id in missing_participants:
            sg_list = index.by_participant[participant_id]

            somalier_paths = {
                info.sg_id: str(all_somalier[info.sg_id]) for info in sg_list if info.sg_id in all_somalier
            }
            if len(somalier_paths) < _MIN_SGS_FOR_IDENTITY_CHECK:
                continue

            jobs = relate.identity_check_jobs(
                participant_id=participant_id,
                outputs=outputs,
                somalier_paths=somalier_paths,
                dataset_name=dataset.name,
                job_attrs={'participant': participant_id},
            )
            all_jobs.extend(jobs)

        return self.make_outputs(dataset, data=outputs, jobs=all_jobs)


@stage.stage(required_stages=[GenerateMissingSomalierFingerprints])
class SomalierPedigreeCheck(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        prefix = dataset.prefix() / 'somalier_checks' / 'pedigree'
        web_prefix = dataset.web_prefix() / 'somalier_checks' / 'pedigree'

        output_prefix = prefix / dataset.name
        web_output_prefix = web_prefix / dataset.name

        return {
            'samples': to_path(f'{output_prefix}.samples.tsv'),
            'pairs': to_path(f'{output_prefix}.pairs.tsv'),
            'expected_ped': prefix / f'{dataset.name}.expected.ped',
            'html': to_path(f'{web_output_prefix}.html'),
            'checks': prefix / f'{dataset.name}-checks.done',
            'output_prefix': str(output_prefix),
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        all_somalier = inputs.as_dict(dataset, GenerateMissingSomalierFingerprints)
        somalier_paths = {sg_id: str(path) for sg_id, path in all_somalier.items()}

        index = get_project_sgs_and_fingerprints(dataset.name)

        # Build PED file content and write to GCS at orchestration time
        ped_content = build_ped_content(dataset.name, index)
        with outputs['expected_ped'].open('w') as f:
            f.write(ped_content)

        jobs = relate.pedigree_check_jobs(
            somalier_paths=somalier_paths,
            outputs=outputs,
            dataset_name=dataset.name,
            label=f'{dataset.name} Somalier',
            job_attrs={},
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)
