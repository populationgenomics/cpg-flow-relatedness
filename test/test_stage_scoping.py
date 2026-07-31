"""
Confirms GenerateMissingSomalierFingerprints is scoped to the whole dataset, not the input cohort.

Previously `_relevant_sg_ids` narrowed the SG set to participants with at least one SG in the input
cohort, so a run over a small cohort only ever extracted fingerprints for that cohort (plus sibling
SGs of the same participants). These tests pin the dataset-wide behaviour: the stage works purely
from the project-wide metamist index, and never consults `dataset.get_sequencing_groups()`.
"""

from typing import NoReturn

import pytest

from rd_qc.stages import GenerateMissingSomalierFingerprints
from rd_qc.utils import SgSomalierInfo, SomalierIndex

# The decorator wraps the class in a factory; __wrapped__ gets the class back so we can build an
# instance without a live Workflow/config (Stage.__init__ needs both).
_STAGE_CLS = GenerateMissingSomalierFingerprints.__wrapped__

_PROJECT = 'my-dataset'

# CPG01 is the only SG in the input cohort. CPG02 shares its participant. CPG03/CPG04 belong to
# participants with no cohort representation at all - these are what the old scoping dropped.
_INDEX_ENTRIES = [
    SgSomalierInfo(sg_id='CPG01', participant_id='PART_A', somalier_path=None),
    SgSomalierInfo(sg_id='CPG02', participant_id='PART_A', somalier_path='gs://bucket/CPG02.cram.somalier'),
    SgSomalierInfo(sg_id='CPG03', participant_id='PART_B', somalier_path=None),
    SgSomalierInfo(sg_id='CPG04', participant_id='PART_C', somalier_path=None),
]

_SOURCE_FILES = {
    'CPG01': 'gs://bucket/CPG01.cram',
    'CPG03': 'gs://bucket/CPG03.cram',
    # CPG04 deliberately absent: no suitable source file for extraction
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

    def _fake_index(project: str) -> SomalierIndex:
        assert project == _PROJECT
        return SomalierIndex(list(_INDEX_ENTRIES))

    def _fake_targets(project: str, sgids: tuple[str, ...]) -> tuple[dict[str, str], dict[str, dict]]:
        assert project == _PROJECT
        calls['requested_sgids'] = sgids
        targets = {sg_id: _SOURCE_FILES[sg_id] for sg_id in sgids if sg_id in _SOURCE_FILES}
        return targets, {sg_id: {} for sg_id in targets}

    def _fake_somalier_jobs(**kwargs: object) -> list:
        calls['job_kwargs'] = kwargs
        return []

    monkeypatch.setattr('rd_qc.stages.get_project_sgs_and_fingerprints', _fake_index)
    monkeypatch.setattr('rd_qc.stages.select_somalier_extract_targets', _fake_targets)
    monkeypatch.setattr('rd_qc.stages.generate_somalier.somalier_jobs', _fake_somalier_jobs)

    instance = object.__new__(_STAGE_CLS)  # bypass Stage.__init__, which needs a live Workflow
    return instance, calls


def test_expected_outputs_cover_every_dataset_sg(stage_instance):
    """Outputs include out-of-cohort participants' SGs, with existing fingerprints reused as-is."""
    instance, _ = stage_instance

    outputs = instance.expected_outputs(_DatasetStub())

    assert {sg_id: str(path) for sg_id, path in outputs.items()} == {
        'CPG01': 'gs://bucket/CPG01.cram.somalier',
        'CPG02': 'gs://bucket/CPG02.cram.somalier',
        'CPG03': 'gs://bucket/CPG03.cram.somalier',
    }


def test_queue_jobs_extracts_for_out_of_cohort_sgs(stage_instance):
    """Extraction is requested for every fingerprint-less SG in the dataset, cohort membership aside."""
    instance, calls = stage_instance

    outputs = instance.queue_jobs(_DatasetStub(), None)

    # CPG02 already has a fingerprint, so it is not re-extracted
    assert calls['requested_sgids'] == ('CPG01', 'CPG03', 'CPG04')

    job_kwargs = calls['job_kwargs']
    assert set(job_kwargs['somalier_targets']) == {'CPG01', 'CPG03'}
    assert set(job_kwargs['somalier_outputs']) == {'CPG01', 'CPG03'}
    assert job_kwargs['project'] == _PROJECT

    # the stage still advertises the full dataset-wide fingerprint set downstream
    assert set(outputs.data) == {'CPG01', 'CPG02', 'CPG03'}


def test_sgs_without_a_source_file_are_dropped(stage_instance):
    """CPG04 has no cram/gvcf to extract from, so it appears in neither outputs nor jobs."""
    instance, calls = stage_instance

    outputs = instance.queue_jobs(_DatasetStub(), None)

    assert 'CPG04' in calls['requested_sgids']
    assert 'CPG04' not in outputs.data
    assert 'CPG04' not in calls['job_kwargs']['somalier_targets']
