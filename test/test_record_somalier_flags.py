"""
Tests for Somalier flag reconciliation against the flags stored in Metamist sequencing-group meta.

Reconciliation overwrites the whole `somalier_flags` list, so these tests focus on what survives
that overwrite — in particular that flag categories absent from the current run are marked resolved
rather than silently dropped, and that the fields derived from the measurement refresh with it.

The Metamist mutation is intercepted rather than sent. Its variable payload is asserted on as the
observable output of reconciliation, since nothing else about a run is visible to a caller.
"""

import pytest

from rd_qc.scripts import record_somalier_flags

FIRST_SEEN = '2026-01-01T00:00:00+00:00'
RESOLVED_EARLIER = '2026-02-02T00:00:00+00:00'
TODAY = '2026-07-30T00:00:00+00:00'

# The SG whose meta is being reconciled. Its pair partner sorts before it, so a key built by
# sorting the pair reads differently from one built by taking the fields in recorded order.
OWNING_SG = 'TST2'
PARTNER_SG = 'TST1'


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
        'sg_id_1': OWNING_SG,
        'sg_id_2': PARTNER_SG,
        'family_external_id': 'FAM1',
        'expected_relationship': 'siblings',
        'inferred_relationship': 'unrelated',
        'relatedness': 0.02,
        'ibs0': 900,
        'ibs2': 100,
        'verdict': 'conflict',
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
        sg={'id': OWNING_SG, 'meta': {'somalier_flags': current_flags}},
        new_flags_by_sg={OWNING_SG: new_flags},
        dataset='test-dataset',
        today=TODAY,
    )


def flags_by_category(written_meta: dict) -> dict[str, dict]:
    """Index the written flags by category — these tests use at most one flag per category."""
    return {flag['category']: flag for flag in written_meta['sgMeta']['somalier_flags']}


# ---------------------------------------------------------------------------
# A category that produced no findings this run
# ---------------------------------------------------------------------------
def test_a_category_absent_from_this_run_is_marked_resolved(written_meta):
    """
    Regression test: the per-category reconcilers used to be skipped when the current run had no new
    flags of that category, so the whole-list overwrite erased those flags with no resolution record.
    """
    reconcile(current_flags=[sex_flag()], new_flags=[relatedness_flag()])

    written = flags_by_category(written_meta)
    assert set(written) == {'sex_inference_mismatch', 'relatedness_mismatch'}
    assert written['sex_inference_mismatch']['resolved'] is True
    assert written['sex_inference_mismatch']['resolution_date'] == TODAY


def test_resolution_preserves_the_first_detected_date(written_meta):
    reconcile(current_flags=[sex_flag()], new_flags=[relatedness_flag()])

    assert flags_by_category(written_meta)['sex_inference_mismatch']['date'] == FIRST_SEEN


def test_this_runs_flags_are_written_unresolved(written_meta):
    reconcile(current_flags=[sex_flag()], new_flags=[relatedness_flag()])

    assert flags_by_category(written_meta)['relatedness_mismatch']['resolved'] is False


def test_already_resolved_flag_of_absent_category_is_kept_untouched(written_meta):
    """An already-resolved flag stays in meta, and its original resolution date is not re-stamped."""
    stored = sex_flag(resolved=True, resolution_date=RESOLVED_EARLIER)

    reconcile(current_flags=[stored], new_flags=[relatedness_flag()])

    kept = flags_by_category(written_meta)['sex_inference_mismatch']
    assert kept['resolved'] is True
    assert kept['resolution_date'] == RESOLVED_EARLIER


# ---------------------------------------------------------------------------
# A finding that is still present
# ---------------------------------------------------------------------------
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


def test_a_retained_flag_takes_the_new_verdict(written_meta):
    # The verdict is derived from the measurement, which drifts between runs while the flag's
    # identity does not. A pair that moves from second-degree into the background stops being a
    # conflict, and a stale 'conflict' would keep it in the section that demands a decision.
    stored = relatedness_flag(verdict='conflict')
    rerun = relatedness_flag(verdict='refinement', relatedness=0.05)

    reconcile(current_flags=[stored], new_flags=[rerun])

    assert flags_by_category(written_meta)['relatedness_mismatch']['verdict'] == 'refinement'


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


# ---------------------------------------------------------------------------
# The SG key stamped onto every flag
# ---------------------------------------------------------------------------
def test_a_pairwise_flag_is_stamped_with_its_sorted_pair_key(written_meta):
    # The report reads this key to dedupe a pair that appears under two families, so it has to be
    # the same string whichever member of the pair is holding the flag.
    reconcile(current_flags=[], new_flags=[relatedness_flag()])

    written = flags_by_category(written_meta)['relatedness_mismatch']

    assert written['sequencing_group_key'] == f'{PARTNER_SG}_{OWNING_SG}'


def test_a_per_sg_flag_is_stamped_with_its_owning_sg(written_meta):
    reconcile(current_flags=[], new_flags=[sex_flag()])

    assert flags_by_category(written_meta)['sex_inference_mismatch']['sequencing_group_key'] == OWNING_SG


def test_a_flag_stored_before_the_key_existed_picks_one_up(written_meta):
    # `relatedness_flag` carries no `sequencing_group_key`, which is exactly the shape of a record
    # written before the field was added. The stored copy must be stamped too, otherwise the
    # report keeps falling back for the lifetime of the flag.
    reconcile(current_flags=[relatedness_flag()], new_flags=[])

    written = flags_by_category(written_meta)['relatedness_mismatch']

    assert written['sequencing_group_key'] == f'{PARTNER_SG}_{OWNING_SG}'


def test_no_stored_or_new_flags_writes_nothing(written_meta):
    """An SG that has never been flagged and is not flagged now should not be mutated at all."""
    reconcile(current_flags=[], new_flags=[])

    assert written_meta == {}
