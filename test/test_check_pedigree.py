"""
Tests for the flag-producing half of the pedigree check.

`produce_flags` is the part that reads somalier's relate outputs and decides what is wrong. It
touches neither Metamist nor Slack, which is what lets the local report harness reuse it, so these
tests drive it directly with synthetic somalier TSVs written to a tmp_path.

The somalier TSV builders live in `fixtures.somalier_tsv`, shared with the consanguinity tests.
"""

import pytest
from fixtures.somalier_tsv import (
    FEMALE,
    IBS0_PARENT_CHILD,
    MALE,
    PROVIDED_SEX_UNKNOWN,
    SEX_INFERENCE_FAILED,
    TRIO_PED,
    pair_row,
    sample_row,
    write_inputs,
)

from rd_qc.scripts import check_pedigree
from rd_qc.scripts.check_pedigree import produce_flags


@pytest.fixture(autouse=True)
def pedigree_messages(monkeypatch) -> list[str]:
    """
    The Slack summary lines one `produce_flags` call appends, isolated per test.

    `check_pedigree._messages` is a module global that `info()` appends to. Swapping in a fresh
    list per test keeps the lines readable here and stops them accumulating across the session.
    """
    lines: list[str] = []
    monkeypatch.setattr(check_pedigree, '_messages', lines)
    return lines


def broken_trio(tmp_path, provided_sex_003: str = 'male') -> dict[str, str]:
    """The trio, measured as three mutually unrelated people. Both parent-child pairs conflict."""
    return write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female'),
            sample_row('TST003', inferred_sex=MALE, provided_sex=provided_sex_003),
        ],
        pair_rows=[
            pair_row('TST001', 'TST002', relatedness=0.01),
            pair_row('TST001', 'TST003', relatedness=0.02),
            pair_row('TST002', 'TST003', relatedness=0.03),
        ],
        ped=TRIO_PED,
    )


def sex_flag_ids(flags: dict) -> set[str]:
    return {sg_id for sg_id, fs in flags.items() for f in fs if f.category == 'sex_inference_mismatch'}


def pedigree_flags(flags: dict) -> list:
    return [f for fs in flags.values() for f in fs if f.category == 'relatedness_mismatch']


# ---------------------------------------------------------------------------
# The happy path
# ---------------------------------------------------------------------------
def test_matching_pedigree_and_sex_produces_no_flags(tmp_path):
    # The measurements confirm the PED: both parent-child pairs at ~0.5 with ibs0 ~0, the two
    # co-parents unrelated, and every sex agreeing.
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female'),
            sample_row('TST003', inferred_sex=MALE, provided_sex='male'),
        ],
        pair_rows=[
            pair_row('TST001', 'TST002', relatedness=0.01),
            pair_row('TST001', 'TST003', relatedness=0.49, ibs0=IBS0_PARENT_CHILD),
            pair_row('TST002', 'TST003', relatedness=0.51, ibs0=IBS0_PARENT_CHILD),
        ],
        ped=TRIO_PED,
    )

    flags, _, _ = produce_flags(**inputs)

    assert flags == {}


# ---------------------------------------------------------------------------
# Sex inference
# ---------------------------------------------------------------------------
def test_sex_mismatch_is_flagged_against_the_owning_sg(tmp_path):
    flags, _, _ = produce_flags(**broken_trio(tmp_path, provided_sex_003='female'))

    sex_flags = [f for f in flags.get('TST003', []) if f.category == 'sex_inference_mismatch']

    assert len(sex_flags) == 1
    assert (sex_flags[0].provided, sex_flags[0].inferred) == ('female', 'male')


def test_matching_sex_is_not_flagged(tmp_path):
    flags, _, _ = produce_flags(**broken_trio(tmp_path, provided_sex_003='male'))

    assert flags, 'the broken trio must still produce its pedigree flags, or this passes vacuously'
    assert sex_flag_ids(flags) == set()


def test_an_unknown_provided_sex_is_not_a_sex_mismatch(tmp_path):
    # A PED that never recorded this person's sex disagrees with the inference without
    # contradicting it. Counting it as a mismatch would flag every dataset with incomplete sex.
    flags, _, _ = produce_flags(**broken_trio(tmp_path, provided_sex_003=PROVIDED_SEX_UNKNOWN))

    assert flags, 'the broken trio must still produce its pedigree flags, or this passes vacuously'
    assert 'TST003' not in sex_flag_ids(flags)


def test_a_failed_sex_inference_is_not_a_sex_mismatch(tmp_path):
    # somalier writes sex=0 when it could not call the sex at all. Nothing was inferred, so there
    # is nothing for the provided sex to contradict.
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female'),
            sample_row('TST003', inferred_sex=SEX_INFERENCE_FAILED, provided_sex='male'),
        ],
        pair_rows=[
            pair_row('TST001', 'TST002', relatedness=0.01),
            pair_row('TST001', 'TST003', relatedness=0.02),
            pair_row('TST002', 'TST003', relatedness=0.03),
        ],
        ped=TRIO_PED,
    )

    flags, _, _ = produce_flags(**inputs)

    assert flags, 'the unrelated measurements must still produce pedigree flags'
    assert 'TST003' not in sex_flag_ids(flags)


# ---------------------------------------------------------------------------
# Pedigree relatedness
# ---------------------------------------------------------------------------
def test_pedigree_mismatch_is_flagged_for_the_disagreeing_pairs(tmp_path):
    flags, _, _ = produce_flags(**broken_trio(tmp_path))

    pedigree = pedigree_flags(flags)
    pairs = {(f.sg_id_1, f.sg_id_2) for f in pedigree}

    # Both parent-child pairs conflict. TST001/TST002 are the co-parents, whom peddy calls
    # 'mom-dad' and which expects unrelated, so the measurement agrees and raises nothing.
    assert pairs == {('TST001', 'TST003'), ('TST002', 'TST003')}
    assert all(f.expected_relationship == 'parent-child' for f in pedigree)
    assert all(f.inferred_relationship == 'unrelated' for f in pedigree)
    assert all(f.verdict == 'conflict' for f in pedigree)
    assert all(f.family_external_id == 'FAM1' for f in pedigree)


def test_pairwise_flags_are_recorded_against_the_first_sg_of_the_pair_only(tmp_path):
    flags, _, _ = produce_flags(**broken_trio(tmp_path))

    owners = {
        (f.sg_id_1, f.sg_id_2): sg_id for sg_id, fs in flags.items() for f in fs if f.category == 'relatedness_mismatch'
    }

    # This is the convention the report's dedup relies on: the lexicographically first SG owns it.
    assert owners == {
        ('TST001', 'TST003'): 'TST001',
        ('TST002', 'TST003'): 'TST002',
    }
    assert 'TST003' not in flags


def test_a_pair_reported_out_of_order_is_normalised_to_the_lexicographic_first(tmp_path):
    # somalier orders pairs.tsv by its own internal sample order, so sample_a can be the
    # lexicographically later of the two. The report's dedup key assumes it never is.
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST003', inferred_sex=MALE, provided_sex='male'),
        ],
        pair_rows=[pair_row('TST003', 'TST001', relatedness=0.02)],
        ped=TRIO_PED,
    )

    flags, _, _ = produce_flags(**inputs)

    flag = pedigree_flags(flags)[0]
    assert (flag.sg_id_1, flag.sg_id_2) == ('TST001', 'TST003')
    assert list(flags) == ['TST001']


def test_same_family_pairs_with_no_recorded_link_are_not_expected_unrelated(tmp_path):
    """
    Two family members with no blood path between them must not be flagged as 'expected unrelated'.

    TST001 and TST002 are both in FAM1 with no parents recorded, so peddy calls them 'unrelated'.
    The genotypes infer full siblings. Before the reframing this was a conflict claiming the
    pedigree said they were unrelated, which it never did.
    """
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            # somalier --infer reconstructs the sibship: both share the same two parents.
            sample_row('TST001', inferred_sex=MALE, provided_sex='male', paternal='TST004', maternal='TST005'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female', paternal='TST004', maternal='TST005'),
            sample_row('TST004', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST005', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[
            pair_row('TST001', 'TST002', relatedness=0.49),
            pair_row('TST001', 'TST004', relatedness=0.48),
            pair_row('TST001', 'TST005', relatedness=0.51),
            pair_row('TST002', 'TST004', relatedness=0.5),
            pair_row('TST002', 'TST005', relatedness=0.49),
            pair_row('TST004', 'TST005', relatedness=0.01),
        ],
        # The expected pedigree knows only that all four are in FAM1, with no links at all.
        ped=(
            'FAM1\tTST001\t0\t0\t1\t1\nFAM1\tTST002\t0\t0\t2\t1\nFAM1\tTST004\t0\t0\t1\t1\nFAM1\tTST005\t0\t0\t2\t1\n'
        ),
    )

    flags, _, _ = produce_flags(**inputs)
    pedigree = pedigree_flags(flags)

    assert pedigree, 'expected the missing links to still be flagged, just not as "unrelated"'
    assert not [f for f in pedigree if f.expected_relationship == 'unrelated']
    assert all(f.expected_relationship == 'related at unknown level' for f in pedigree)


def test_cross_family_pairs_keep_expected_unrelated(tmp_path):
    """
    Two people in different families really are expected to be unrelated.

    So the reframing must not touch them: a related inference across a family boundary is the
    cross-family swap case, and it has to stay a conflict. Here somalier infers TST001 and TST002
    as full siblings while the expected pedigree has them in FAM1 and FAM2 respectively.
    """
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male', paternal='TST004', maternal='TST005'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female', paternal='TST004', maternal='TST005'),
            sample_row('TST004', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST005', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[pair_row('TST001', 'TST002', relatedness=0.47)],
        ped=(
            'FAM1\tTST001\t0\t0\t1\t1\nFAM2\tTST002\t0\t0\t2\t1\nFAM1\tTST004\t0\t0\t1\t1\nFAM1\tTST005\t0\t0\t2\t1\n'
        ),
    )

    flags, _, _ = produce_flags(**inputs)
    pedigree = pedigree_flags(flags)

    assert [(f.expected_relationship, f.inferred_relationship) for f in pedigree] == [('unrelated', 'siblings')]


# ---------------------------------------------------------------------------
# Zero-depth samples
# ---------------------------------------------------------------------------
def unusable_third_sample(tmp_path) -> dict[str, str]:
    """The trio with TST003 ungenotypable: gt_depth_mean 0 means its calls are meaningless."""
    return write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female'),
            sample_row('TST003', inferred_sex=MALE, provided_sex='female', gt_depth_mean=0.0),
        ],
        pair_rows=[
            pair_row('TST001', 'TST002', relatedness=0.01),
            pair_row('TST001', 'TST003', relatedness=0.02),
            pair_row('TST002', 'TST003', relatedness=0.03),
        ],
        ped=TRIO_PED,
    )


def test_a_zero_depth_sample_gets_no_sex_flag(tmp_path):
    # TST003's provided female against an inferred male is a real disagreement, but the inference
    # came off meaningless genotypes, so reporting it would be a false positive.
    flags, _, _ = produce_flags(**unusable_third_sample(tmp_path))

    assert sex_flag_ids(flags) == set()


def test_no_pair_involving_a_zero_depth_sample_is_flagged(tmp_path):
    flags, _, _ = produce_flags(**unusable_third_sample(tmp_path))

    flagged_pairs = {(f.sg_id_1, f.sg_id_2) for f in pedigree_flags(flags)}

    # Both TST003 pairs measure unrelated against a recorded parent-child link, so they would be
    # flagged on their numbers alone; only the depth exclusion keeps them out.
    assert flagged_pairs == set()


# ---------------------------------------------------------------------------
# The Slack summary's reporting buckets
# ---------------------------------------------------------------------------
def test_an_identical_pair_is_reported_in_the_most_serious_bucket(tmp_path, pedigree_messages):
    # One genome recorded as two people is a duplicate or a swap, never a pedigree omission.
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST003', inferred_sex=MALE, provided_sex='male'),
        ],
        pair_rows=[pair_row('TST001', 'TST003', relatedness=0.98, ibs0=IBS0_PARENT_CHILD)],
        ped=TRIO_PED,
    )

    produce_flags(**inputs)

    assert '1 sample pair(s) inferred as identical' in '\n'.join(pedigree_messages)


def test_a_relationship_the_measurement_lost_is_reported_as_less_related(tmp_path, pedigree_messages):
    produce_flags(**broken_trio(tmp_path))

    # Both parent-child pairs measured unrelated: the sample-swap shape.
    assert '2 sample pair(s) that are recorded as related, but measured as less related' in '\n'.join(pedigree_messages)


def test_a_relationship_the_measurement_added_is_reported_as_more_related(tmp_path, pedigree_messages):
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[pair_row('TST001', 'TST002', relatedness=0.47)],
        ped='FAM1\tTST001\t0\t0\t1\t1\nFAM2\tTST002\t0\t0\t2\t1\n',
    )

    produce_flags(**inputs)

    assert '1 sample pair(s) that are recorded as unrelated, but measured as related' in '\n'.join(pedigree_messages)


def test_refinements_are_counted_rather_than_listed(tmp_path, pedigree_messages):
    # These are routinely the bulk of the flags, and listing them drowns the buckets that need a
    # decision, so the summary states the count and points at the report.
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[pair_row('TST001', 'TST002', relatedness=0.47)],
        ped='FAM1\tTST001\t0\t0\t1\t1\nFAM1\tTST002\t0\t0\t2\t1\n',
    )

    produce_flags(**inputs)

    text = '\n'.join(pedigree_messages)
    assert '1 pair(s) measured as related where the pedigree records no relationship' in text
    # The pair itself is not enumerated under a heading of its own.
    assert 'sample pair(s) that are recorded as' not in text


def test_a_clean_pedigree_reports_the_all_clear(tmp_path, pedigree_messages):
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST003', inferred_sex=MALE, provided_sex='male'),
        ],
        pair_rows=[pair_row('TST001', 'TST003', relatedness=0.49, ibs0=IBS0_PARENT_CHILD)],
        ped=TRIO_PED,
    )

    produce_flags(**inputs)

    assert 'Measured relatedness matches the pedigree for every pair' in '\n'.join(pedigree_messages)
