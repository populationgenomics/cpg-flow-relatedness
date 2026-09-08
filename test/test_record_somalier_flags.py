"""
Tests for Somalier flag reconciliation against the flags stored in Metamist sequencing-group meta.

Reconciliation overwrites the whole `somalier_flags` list, so these tests focus on what survives
that overwrite — in particular that flag categories absent from the current run are marked resolved
rather than silently dropped.
"""

import pytest

from rd_qc.scripts import record_somalier_flags

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

    monkeypatch.setattr(record_somalier_flags, 'query', fake_query)
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
    assert len(written) == 2  # noqa: PLR2004

    by_inferred = {flag['inferred_relationship']: flag for flag in written}
    assert by_inferred['unrelated']['resolved'] is True
    assert by_inferred['unrelated']['resolution_date'] == TODAY
    assert by_inferred['parent-child']['resolved'] is False


def test_no_stored_or_new_flags_writes_nothing(written_meta):
    """An SG that has never been flagged and is not flagged now should not be mutated at all."""
    reconcile(current_flags=[], new_flags=[])

    assert written_meta == {}
