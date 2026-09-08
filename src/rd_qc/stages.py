"""Stages for the rd_qc somalier QC workflow."""

from datetime import datetime

from rd_qc.jobs import generate_somalier, relate, somalier_flags_report
from rd_qc.utils import (
    build_ped_content,
    find_sgids_without_somalier,
    get_project_sgs_and_fingerprints,
    select_somalier_extract_targets,
    sg_ids_tag,
)

from cpg_flow import stage, targets
from cpg_utils import Path, config, to_path
from cpg_utils.existence_checks import exists

_MIN_SGS_FOR_IDENTITY_CHECK = 2

def convert_to_web_url(path: Path, dataset: targets.Dataset) -> str:
    """Convert a Path to a web URL, if the dataset has a web URL."""
    if base_url := dataset.web_url():
        return str(path).replace(str(dataset.web_prefix()), base_url)
    return str(path)

@stage.stage()
class GenerateMissingSomalierFingerprints(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        dataset_name = config.dataset_for_access_level(dataset.name)
        index = get_project_sgs_and_fingerprints(dataset_name)
        outputs: dict[str, Path] = {}

        for info in index.by_sg.values():
            if info.somalier_path is not None:
                outputs[info.sg_id] = to_path(info.somalier_path)

        missing_sgids = find_sgids_without_somalier(index)
        if missing_sgids:
            extract_targets, _ = select_somalier_extract_targets(
                dataset_name,
                tuple(sorted(missing_sgids)),
            )
            for sg_id, source_file in extract_targets.items():
                outputs[sg_id] = to_path(f'{source_file}.somalier')

        return outputs

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:  # noqa: ARG002
        outputs = self.expected_outputs(dataset)

        dataset_name = config.dataset_for_access_level(dataset.name)
        index = get_project_sgs_and_fingerprints(dataset_name)
        missing_sgids = find_sgids_without_somalier(index)

        if not missing_sgids:
            return self.make_outputs(dataset, data=outputs)

        extract_targets, sg_id_map = select_somalier_extract_targets(
            dataset_name,
            tuple(sorted(missing_sgids)),
        )

        # Only pass the missing SGs' output paths to the job builder
        somalier_outputs = {sg_id: outputs[sg_id] for sg_id in extract_targets}

        jobs = generate_somalier.somalier_jobs(
            dataset_name=dataset_name,
            sg_id_map=sg_id_map,
            somalier_targets=extract_targets,
            somalier_outputs=somalier_outputs,
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage.stage(required_stages=[GenerateMissingSomalierFingerprints])
class SomalierSelfCheck(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        dataset_name = config.dataset_for_access_level(dataset.name)
        index = get_project_sgs_and_fingerprints(dataset_name)
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
            outputs[f'{participant_id}_json'] = to_path(str(prefix) + '.checks.json')

        return outputs

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        if not outputs:
            return self.make_outputs(dataset, data=outputs)

        dataset_name = config.dataset_for_access_level(dataset.name)
        all_somalier = inputs.as_dict(dataset, GenerateMissingSomalierFingerprints)
        index = get_project_sgs_and_fingerprints(dataset_name)

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

            jobs = relate.self_relatedness_jobs(
                dataset_name=dataset_name,
                participant_id=participant_id,
                somalier_paths=somalier_paths,
                outputs=outputs,
                job_attrs={'participant': participant_id},
            )
            all_jobs.extend(jobs)

        return self.make_outputs(dataset, data=outputs, jobs=all_jobs)


@stage.stage(required_stages=[GenerateMissingSomalierFingerprints, SomalierSelfCheck])
class SomalierPedigreeCheck(stage.DatasetStage):
    """
    This stage runs somalier relate on all SGs in the dataset using the expected pedigree. It then compiles these
    results alongside as the self-relatedness checks from the previous stage and registers any flags in metamist.

    Flags are determined by comparing the relatedness results to the relationships defined in the expected pedigree.

    While the previous stages are run on all SGs in the dataset, this stage is scoped to only those SGs that meet
    the sequencing type & technology requirements as defined in the config. This prevents the results from becoming
    too difficult to interpret for datasets where participants have SGs with multiple different sequencing types
    and technologies.
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        """
        Expected outputs for the somalier pedigree check stage.

        Files are written to paths namespaced by the AR GUID to avoid collisions between runs of the workflow.

        The final HTML report is written to both a path with and without namespacing, so that it is updated
        with each run, whilst also preserving the previous run's report for reference.
        """
        ar_guid: str = config.config_retrieve(['workflow', 'ar-guid'])
        prefix = dataset.prefix() / 'somalier_checks' / 'pedigree' / ar_guid
        output_prefix = prefix / dataset.name

        base_web_prefix = dataset.web_prefix() / 'somalier_checks' / 'pedigree'
        base_web_output_prefix = base_web_prefix / dataset.name

        web_prefix = base_web_prefix / ar_guid
        web_output_prefix = web_prefix / dataset.name

        return {
            'samples': to_path(f'{output_prefix}.samples.tsv'),
            'pairs': to_path(f'{output_prefix}.pairs.tsv'),
            'expected_ped': prefix / f'{dataset.name}.expected.ped',
            'html': to_path(f'{web_output_prefix}.html'),
            'base_html_url': to_path(f'{base_web_output_prefix}.html'),
            'checks': prefix / f'{dataset.name}-checks.done',
            'json': to_path(f'{output_prefix}.checks.json'),
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        # filter_sgs=True ensures that only SGs meeting the sequencing type & technology requirements are included
        index = get_project_sgs_and_fingerprints(dataset.name, filter_sgs=True)
        somalier_paths = {sg_id: info.somalier_path for sg_id, info in index.by_sg.items()}

        somalier_self_relatedness_json_paths = [
            path for path in inputs.as_dict(dataset, SomalierSelfCheck).values() if str(path).endswith('.json')
        ]

        # Build PED file content and write to GCS at orchestration time
        ped_content = build_ped_content(dataset.name, index)
        with outputs['expected_ped'].open('w') as f:
            f.write(ped_content)

        jobs = relate.pedigree_check_jobs(
            somalier_paths=somalier_paths,
            somalier_self_relatedness_json_paths=somalier_self_relatedness_json_paths,
            outputs=outputs,
            tmp_prefix=dataset.tmp_prefix() / 'somalier_checks' / 'pedigree',
            dataset_name=dataset.name,
            label=f'{dataset.name} Somalier',
            job_attrs={},
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)

@stage.stage(
    required_stages=[SomalierPedigreeCheck],
    forced=True,
)
class GenerateSomalierFlagsReport(stage.DatasetStage):
    """
    Queries Metamist for all Somalier flags across the dataset's sequencing groups
    and generates a summary HTML report saved to both a static URL and a
    timestamped URL in the dataset's web bucket.
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        timestamp = datetime.now().astimezone().strftime('%Y-%m-%d_%H%M%S')
        return {
            'timestamped': dataset.web_prefix() / 'somalier_flags' / timestamp / 'somalier_flags_report.html',
            'html': dataset.web_prefix() / 'somalier_flags' / 'somalier_flags_report.html',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        out_html_url = convert_to_web_url(outputs['html'], dataset)
        somalier_report_url = convert_to_web_url(
            inputs.as_path_by_target(SomalierPedigreeCheck, 'base_html_url')[dataset.name], dataset)

        jobs = somalier_flags_report.somalier_flags_report_job(
            dataset=dataset.name,
            outputs=outputs,
            out_html_url=out_html_url,
            somalier_report_url=somalier_report_url,
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=outputs, jobs=jobs)
