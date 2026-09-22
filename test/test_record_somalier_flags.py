"""
Tests for Somalier flag reconciliation against the flags stored in Metamist sequencing-group meta.

Reconciliation overwrites the whole `somalier_flags` list, so these tests focus on what survives
that overwrite — in particular that flag categories absent from the current run are marked resolved
rather than silently dropped.
"""

import pytest

from rd_qc import flag_store
from rd_qc.scripts import record_somalier_flags
from rd_qc.utils import SomalierRelatednessFlag

FIRST_SEEN = '2026-01-01T00:00:00+00:00'
RESOLVED_EARLIER = '2026-02-02T00:00:00+00:00'
TODAY = '2026-07-30T00:00:00+00:00'


def sex_flag(**overrides: object) -> dict:
    """A sex_inference_mismatch flag as stored in / written to SG meta."""
    return {
        'category': 'sex_inference_mismatch',
        'date': FIRST_SEEN,
        'ar_guid': 'guid-old',
        'resolved': False,
        'resolution_date': None,
        'provided': 'M',
        'inferred': 'F',
        'mean_depth': 30.0,
        'x_het_ratio': 0.4,
        'x_depth_ratio': 2.0,
        'y_depth_ratio': 0.0,
        'x_sites': 500,
        'p_middling_ab': 0.01,
    } | overrides


def self_relatedness_flag(**overrides: object) -> dict:
    """A self_relatedness_mismatch flag as stored in / written to SG meta."""
    return {
        'category': 'self_relatedness_mismatch',
        'date': FIRST_SEEN,
        'ar_guid': 'guid-old',
        'resolved': False,
        'resolution_date': None,
        'sg_id_1': 'CPG1',
        'sg_id_2': 'CPG1B',
        'participant_external_id': 'PART1',
        'threshold': 0.4,
        'relatedness': 0.9,
        'ibs0': 5,
        'ibs2': 950,
    } | overrides


def relatedness_flag(**overrides: object) -> dict:
    """A relatedness_mismatch flag as stored in / written to SG meta."""
    return {
        'category': 'relatedness_mismatch',
        'date': FIRST_SEEN,
        'ar_guid': 'guid-old',
        'resolved': False,
        'resolution_date': None,
        'sg_id_1': 'CPG1',
        'sg_id_2': 'CPG2',
        'family_external_id': 'FAM1',
        'expected_relationship': 'siblings',
        'inferred_relationship': 'unrelated',
        'relatedness': 0.02,
        'ibs0': 900,
        'ibs2': 100,
    } | overrides


@pytest.fixture
def written_meta(monkeypatch):
    """
    Intercept the Metamist mutation and expose the meta payload it would have written.

    Returns a dict that stays empty if no mutation was attempted.
    """
    captured: dict = {}

    def fake_query(_query, variables=None) -> dict:
        captured.update(variables or {})
        return {}

    monkeypatch.setattr(flag_store, 'query', fake_query)
    return captured


def reconcile(current_flags: list[dict], new_flags: list[dict]) -> None:
    """Reconcile one SG's flags, given what is stored and what this run produced."""
    record_somalier_flags.reconcile_sg_somalier_flags(
        sg={'id': 'CPG1', 'meta': {'somalier_flags': current_flags}},
        new_flags_by_sg={'CPG1': new_flags},
        dataset='test-dataset',
        today=TODAY,
    )


def flags_by_category(written_meta: dict) -> dict[str, dict]:
    """Index the written flags by category — these tests use at most one flag per category."""
    return {flag['category']: flag for flag in written_meta['sgMeta']['somalier_flags']}


def test_absent_category_is_resolved_not_dropped(written_meta):
    """
    A stored flag whose category produced no findings this run must be retained and marked resolved.

    Regression test: the per-category reconcilers used to be skipped when the current run had no new
    flags of that category, so the whole-list overwrite erased those flags with no resolution record.
    """
    reconcile(current_flags=[sex_flag()], new_flags=[relatedness_flag()])

    written = flags_by_category(written_meta)
    assert set(written) == {'sex_inference_mismatch', 'relatedness_mismatch'}

    resolved_sex = written['sex_inference_mismatch']
    assert resolved_sex['resolved'] is True
    assert resolved_sex['resolution_date'] == TODAY
    assert resolved_sex['date'] == FIRST_SEEN, 'first-detected date must survive resolution'

    assert written['relatedness_mismatch']['resolved'] is False


def test_already_resolved_flag_of_absent_category_is_kept_untouched(written_meta):
    """An already-resolved flag stays in meta, and its original resolution date is not re-stamped."""
    stored = sex_flag(resolved=True, resolution_date=RESOLVED_EARLIER)

    reconcile(current_flags=[stored], new_flags=[relatedness_flag()])

    kept = flags_by_category(written_meta)['sex_inference_mismatch']
    assert kept['resolved'] is True
    assert kept['resolution_date'] == RESOLVED_EARLIER


def test_recurring_flag_is_retained_with_refreshed_measurements(written_meta):
    """
    A flag with unchanged identity is retained, keeping its original date but taking the new metrics.

    `relatedness`/`ibs0`/`ibs2` drift between relate runs, so they are excluded from flag identity.
    """
    recurrence = relatedness_flag(date=TODAY, ar_guid='guid-new', relatedness=0.05, ibs0=850, ibs2=120)

    reconcile(current_flags=[relatedness_flag()], new_flags=[recurrence])

    written = flags_by_category(written_meta)
    assert set(written) == {'relatedness_mismatch'}

    retained = written['relatedness_mismatch']
    assert retained['resolved'] is False
    assert retained['resolution_date'] is None
    assert retained['date'] == FIRST_SEEN, 'a recurring issue keeps its first-detected date'
    assert (retained['relatedness'], retained['ibs0'], retained['ibs2']) == (0.05, 850, 120)


def test_changed_identity_resolves_the_old_flag_and_adds_a_new_one(written_meta):
    """A different inferred relationship is a different finding, so the previous one is resolved."""
    reinferred = relatedness_flag(inferred_relationship='parent-child')

    reconcile(current_flags=[relatedness_flag()], new_flags=[reinferred])

    written = written_meta['sgMeta']['somalier_flags']
    assert len(written) == 2

    by_inferred = {flag['inferred_relationship']: flag for flag in written}
    assert by_inferred['unrelated']['resolved'] is True
    assert by_inferred['unrelated']['resolution_date'] == TODAY
    assert by_inferred['parent-child']['resolved'] is False


def test_no_stored_or_new_flags_writes_nothing(written_meta):
    """An SG that has never been flagged and is not flagged now should not be mutated at all."""
    reconcile(current_flags=[], new_flags=[])

    assert written_meta == {}


def test_manual_resolution_fields_default_to_absent():
    """A flag nobody has reviewed carries the fields, unset, so every record has the same shape."""
    flag = SomalierRelatednessFlag(
        category='relatedness_mismatch',
        sg_id_1='CPG1',
        sg_id_2='CPG2',
        family_external_id='FAM1',
        expected_relationship='siblings',
        inferred_relationship='unrelated',
        relatedness=0.02,
        ibs0=900,
        ibs2=100,
    )

    assert flag.manually_resolved is False
    assert flag.manual_resolution_reason is None
    assert flag.manual_resolution_by is None


def test_a_flag_stored_before_the_manual_fields_existed_still_deserialises():
    """The fixture dict has none of the new keys, which is what Metamist holds for older flags."""
    flag = SomalierRelatednessFlag(**relatedness_flag())

    assert flag.manually_resolved is False
    assert flag.manual_resolution_reason is None
    assert flag.manual_resolution_by is None


def curator_resolved(flag: dict, **overrides: object) -> dict:
    """`flag` as the resolve CLI leaves it: resolved, with the reviewer's reason attached."""
    return (
        flag
        | {
            'resolved': True,
            'resolution_date': RESOLVED_EARLIER,
            'manually_resolved': True,
            'manual_resolution_reason': 'pedigree known wrong',
            'manual_resolution_by': 'ef',
        }
        | overrides
    )


def test_manually_resolved_flag_stays_resolved_when_the_finding_recurs(written_meta):
    """
    The whole point of the feature: a run that measures the same thing again must not reopen it.

    Without the manual branch this flag falls through compare_* into the overwrite branch, which
    sets resolved=False and leaves the manual fields on an active flag.
    """
    held = curator_resolved(relatedness_flag())
    recurrence = relatedness_flag(date=TODAY, relatedness=0.05, ibs0=850, ibs2=120)

    reconcile(current_flags=[held], new_flags=[recurrence])

    written = flags_by_category(written_meta)['relatedness_mismatch']
    assert written['resolved'] is True
    assert written['manually_resolved'] is True
    assert written['manual_resolution_by'] == 'ef'
    assert written['manual_resolution_reason'] == 'pedigree known wrong'
    assert written['resolution_date'] == RESOLVED_EARLIER, 'the reviewer resolved it, not this run'
    assert written['date'] == FIRST_SEEN, 'a held issue keeps its first-detected date'
    assert (written['relatedness'], written['ibs0'], written['ibs2']) == (0.05, 850, 120)


def test_a_changed_finding_is_not_suppressed_by_a_manual_resolution(written_meta):
    """
    Binding is strict: the marker is on one record, so a different measurement surfaces unheld.

    This is the safety property. A pair accepted as parent-child must not stay quiet when the
    genotypes start saying unrelated.
    """
    held = curator_resolved(relatedness_flag())
    reinferred = relatedness_flag(inferred_relationship='parent-child')

    reconcile(current_flags=[held], new_flags=[reinferred])

    written = written_meta['sgMeta']['somalier_flags']
    assert len(written) == 2

    by_inferred = {flag['inferred_relationship']: flag for flag in written}
    assert by_inferred['unrelated']['manually_resolved'] is True
    assert by_inferred['parent-child']['resolved'] is False
    assert by_inferred['parent-child']['manually_resolved'] is False


def test_a_manually_resolved_flag_stays_held_when_the_finding_disappears(written_meta):
    """An accepted finding that later goes away keeps its marker and stays out of the report."""
    held = curator_resolved(relatedness_flag())

    reconcile(current_flags=[held], new_flags=[sex_flag()])

    written = flags_by_category(written_meta)['relatedness_mismatch']
    assert written['manually_resolved'] is True
    assert written['resolution_date'] == RESOLVED_EARLIER, 'the resolution date is not re-stamped'


@pytest.mark.parametrize(
    ('category', 'flag_factory', 'recurrence_overrides', 'measured_field', 'refreshed_value'),
    [
        ('sex_inference_mismatch', sex_flag, {'mean_depth': 29.0}, 'mean_depth', 29.0),
        ('self_relatedness_mismatch', self_relatedness_flag, {'relatedness': 0.5}, 'relatedness', 0.5),
        ('relatedness_mismatch', relatedness_flag, {'relatedness': 0.05}, 'relatedness', 0.05),
    ],
)
def test_manual_resolution_is_held_for_every_category(
    written_meta, category, flag_factory, recurrence_overrides, measured_field, refreshed_value
):
    """
    The branch is duplicated across three reconcilers, so all three need their own coverage.

    Each case recurs the SAME identity with only a measured field changed, so this exercises the
    held branch specifically rather than the added-new-flag path (which a changed identity would
    trigger instead, passing for the wrong reason).
    """
    held = curator_resolved(flag_factory())

    reconcile(current_flags=[held], new_flags=[flag_factory(**recurrence_overrides)])

    written = flags_by_category(written_meta)[category]
    assert written['resolved'] is True
    assert written['manually_resolved'] is True
    assert written[measured_field] == refreshed_value, 'measured values still refresh while held'


@pytest.mark.parametrize(
    ('category', 'factory', 'refreshed', 'field', 'expected'),
    [
        ('sex_inference_mismatch', sex_flag, {'mean_depth': 29.0}, 'mean_depth', 29.0),
        ('self_relatedness_mismatch', self_relatedness_flag, {'relatedness': 0.5}, 'relatedness', 0.5),
        ('relatedness_mismatch', relatedness_flag, {'relatedness': 0.05}, 'relatedness', 0.05),
    ],
)
def test_a_recurring_flag_is_retained_for_every_category(written_meta, category, factory, refreshed, field, expected):
    """
    The retained branch, pinned per category.

    All three reconcilers route this branch through one `refresh_measured_values` call with their
    own field tuple, so passing the wrong tuple would silently stop refreshing measurements. Only
    the relatedness category covered this before, which is how a wrong tuple could have shipped.
    """
    reconcile(current_flags=[factory()], new_flags=[factory(date=TODAY, **refreshed)])

    written = flags_by_category(written_meta)[category]
    assert written['resolved'] is False
    assert written['resolution_date'] is None
    assert written['date'] == FIRST_SEEN, 'a recurring issue keeps its first-detected date'
    assert written[field] == expected, "this run's measurement is taken"
