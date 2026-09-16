"""
Confirms GenerateMissingSomalierFingerprints is scoped to the whole dataset, not the input cohort.

Previously `_relevant_sg_ids` narrowed the SG set to participants with at least one SG in the input
cohort, so a run over a small cohort only ever extracted fingerprints for that cohort (plus sibling
SGs of the same participants). These tests pin the dataset-wide behaviour: the stage works purely
from the project-wide metamist index, and never consults `dataset.get_sequencing_groups()`.

Access-level resolution is supplied by the fixture rather than read from ambient config, so the
suite's result does not depend on whether the developer has CPG_CONFIG_PATH set.
"""

from typing import NoReturn

import pytest

from rd_qc.stages import GenerateMissingSomalierFingerprints
from rd_qc.utils import SgSomalierInfo, SomalierIndex

from cpg_utils import config

# The decorator wraps the class in a factory; __wrapped__ gets the class back so we can build an
# instance without a live Workflow/config (Stage.__init__ needs both).
_STAGE_CLS = GenerateMissingSomalierFingerprints.__wrapped__

# The dataset as cpg-flow names it, and the project name access-level resolution turns it into.
# They differ so a stage that queries metamist with the unresolved name is caught.
_PROJECT = 'my-dataset'
_RESOLVED_PROJECT = 'my-dataset-test'

# TST01 is the only SG in the input cohort. TST02 shares its participant. TST03/TST04 belong to
# participants with no cohort representation at all - these are what the old scoping dropped.
_INDEX_ENTRIES = [
    SgSomalierInfo(sg_id='TST01', participant_id=1, participant_external_id='PART_A', somalier_path=None),
    SgSomalierInfo(
        sg_id='TST02',
        participant_id=1,
        participant_external_id='PART_A',
        somalier_path='gs://bucket/TST02.cram.somalier',
    ),
    SgSomalierInfo(sg_id='TST03', participant_id=2, participant_external_id='PART_B', somalier_path=None),
    SgSomalierInfo(sg_id='TST04', participant_id=3, participant_external_id='PART_C', somalier_path=None),
]

_SOURCE_FILES = {
    'TST01': 'gs://bucket/TST01.cram',
    'TST03': 'gs://bucket/TST03.cram',
    # TST04 deliberately absent: no suitable source file for extraction
}


class _DatasetStub:
    """Stands in for cpg_flow.targets.Dataset, which only ever holds the input cohorts' SGs."""

    name = _PROJECT

    def get_sequencing_groups(self) -> NoReturn:
        raise AssertionError('stage must not scope itself to the input cohort SGs')


@pytest.fixture
def stage_instance(monkeypatch):
    """A stage instance plus a record of what it asked metamist and the job builder for."""
    calls: dict[str, object] = {}

    def _resolve_access_level(dataset_name: str) -> str:
        """Stands in for ambient-config access-level resolution. Raises rather than guessing."""
        if dataset_name != _PROJECT:
            raise AssertionError(f'unexpected dataset name {dataset_name!r}')
        return _RESOLVED_PROJECT

    def _fake_index(project: str) -> SomalierIndex:
        calls['index_project'] = project
        return SomalierIndex(list(_INDEX_ENTRIES))

    def _fake_targets(project: str, sgids: tuple[str, ...]) -> tuple[dict[str, str], dict[str, dict]]:
        calls['targets_project'] = project
        calls['requested_sgids'] = sgids
        targets = {sg_id: _SOURCE_FILES[sg_id] for sg_id in sgids if sg_id in _SOURCE_FILES}
        return targets, {sg_id: {} for sg_id in targets}

    def _fake_somalier_jobs(**kwargs: object) -> list:
        calls['job_kwargs'] = kwargs
        return []

    monkeypatch.setattr(config, 'dataset_for_access_level', _resolve_access_level)
    monkeypatch.setattr('rd_qc.stages.get_project_sgs_and_fingerprints', _fake_index)
    monkeypatch.setattr('rd_qc.stages.select_somalier_extract_targets', _fake_targets)
    monkeypatch.setattr('rd_qc.stages.generate_somalier.somalier_jobs', _fake_somalier_jobs)

    instance = object.__new__(_STAGE_CLS)  # bypass Stage.__init__, which needs a live Workflow
    return instance, calls


# ---------------------------------------------------------------------------
# expected_outputs
# ---------------------------------------------------------------------------
def test_expected_outputs_cover_every_dataset_sg(stage_instance):
    """Outputs include out-of-cohort participants' SGs, not just the input cohort's."""
    instance, _ = stage_instance

    outputs = instance.expected_outputs(_DatasetStub())

    # TST04 has no source file to extract from, so it cannot be promised.
    assert set(outputs) == {'TST01', 'TST02', 'TST03'}


def test_an_existing_fingerprint_is_reused_rather_than_re_derived(stage_instance):
    instance, _ = stage_instance

    outputs = instance.expected_outputs(_DatasetStub())

    assert str(outputs['TST02']) == 'gs://bucket/TST02.cram.somalier'


def test_a_missing_fingerprint_is_named_after_its_source_file(stage_instance):
    instance, _ = stage_instance

    outputs = instance.expected_outputs(_DatasetStub())

    assert str(outputs['TST03']) == 'gs://bucket/TST03.cram.somalier'


# ---------------------------------------------------------------------------
# queue_jobs
# ---------------------------------------------------------------------------
def test_extraction_is_requested_for_every_fingerprintless_sg(stage_instance):
    """Cohort membership aside: TST03 and TST04 have no SG in the input cohort at all."""
    instance, calls = stage_instance

    instance.queue_jobs(_DatasetStub(), None)

    # TST02 already has a fingerprint, so it is not re-extracted.
    assert set(calls['requested_sgids']) == {'TST01', 'TST03', 'TST04'}


def test_extract_targets_are_requested_in_a_stable_order(stage_instance):
    # select_somalier_extract_targets is @cache'd on its arguments, so an unordered set would
    # produce a different cache key per run and re-query metamist every time.
    instance, calls = stage_instance

    instance.queue_jobs(_DatasetStub(), None)

    assert list(calls['requested_sgids']) == sorted(calls['requested_sgids'])


def test_only_the_missing_fingerprints_are_handed_to_the_job_builder(stage_instance):
    instance, calls = stage_instance

    instance.queue_jobs(_DatasetStub(), None)

    job_kwargs = calls['job_kwargs']
    assert set(job_kwargs['somalier_targets']) == {'TST01', 'TST03'}
    assert set(job_kwargs['somalier_outputs']) == {'TST01', 'TST03'}


def test_the_stage_advertises_the_whole_dataset_downstream(stage_instance):
    # Downstream stages relate every SG in the dataset, not only the ones extracted this run.
    instance, _ = stage_instance

    outputs = instance.queue_jobs(_DatasetStub(), None)

    assert set(outputs.data) == {'TST01', 'TST02', 'TST03'}


def test_sgs_without_a_source_file_are_dropped(stage_instance):
    """TST04 has no cram/gvcf to extract from, so it appears in neither outputs nor jobs."""
    instance, calls = stage_instance

    outputs = instance.queue_jobs(_DatasetStub(), None)

    assert 'TST04' in calls['requested_sgids'], 'it must be asked for before it can be dropped'
    assert 'TST04' not in outputs.data
    assert 'TST04' not in calls['job_kwargs']['somalier_targets']


# ---------------------------------------------------------------------------
# Access-level resolution
# ---------------------------------------------------------------------------
def test_metamist_is_queried_with_the_access_level_resolved_project(stage_instance):
    # Querying the unresolved name reads the production project from a test run.
    instance, calls = stage_instance

    instance.queue_jobs(_DatasetStub(), None)

    assert calls['index_project'] == _RESOLVED_PROJECT
    assert calls['targets_project'] == _RESOLVED_PROJECT


def test_the_job_builder_is_given_the_resolved_project(stage_instance):
    instance, calls = stage_instance

    instance.queue_jobs(_DatasetStub(), None)

    assert calls['job_kwargs']['dataset_name'] == _RESOLVED_PROJECT
