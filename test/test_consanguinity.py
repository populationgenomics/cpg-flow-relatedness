"""
Tests for excusing a co-parent relatedness flag with a recorded consanguineous union.

Metamist records consanguinity against the child, in the free-text participant phenotypes, so the
check has to walk from a flagged mom-dad pair to their children. Two halves are tested here: the
extraction of the phenotype out of a Metamist response, and the pedigree walk inside produce_flags.

Absence is deliberately treated as "not recorded", so a dataset that has never supplied the field
keeps every co-parent flag rather than having them all quietly excused.
"""

from pathlib import Path

from test_check_pedigree import (
    FEMALE,
    MALE,
    sample_row,
    write_inputs,
)

from rd_qc.scripts.check_pedigree import produce_flags
from rd_qc.utils import consanguineous_sg_ids

# Second-degree relatedness: above SECOND_DEGREE_MIN_RELATEDNESS, below first-degree.
CONSANGUINEOUS_KIN = 0.21
CONSANGUINEOUS_IBS0 = 800


def sg_response(sg_id: str, phenotypes: dict | None) -> dict:
    """One sequencingGroups entry in the shape SG_QUERY returns."""
    return {'id': sg_id, 'sample': {'participant': {'phenotypes': phenotypes}}}


# ---------------------------------------------------------------------------
# Reading the phenotype out of Metamist's response
# ---------------------------------------------------------------------------
def test_the_recorded_string_one_counts():
    assert consanguineous_sg_ids([sg_response('CPG001', {'Consanguinity': '1'})]) == {'CPG001'}


def test_zero_and_absent_are_both_treated_as_not_recorded():
    sgs = [
        sg_response('CPG001', {'Consanguinity': '0'}),
        sg_response('CPG002', {'Birth Year': '1984'}),
        sg_response('CPG003', {}),
        sg_response('CPG004', None),
    ]

    assert consanguineous_sg_ids(sgs) == set()


def test_a_missing_participant_or_sample_does_not_raise():
    assert consanguineous_sg_ids([{'id': 'CPG001'}, {'id': 'CPG002', 'sample': None}]) == set()


def test_the_key_is_matched_regardless_of_case():
    # 'Consanguinity' is the house spelling, but the phenotypes dict is free text.
    assert consanguineous_sg_ids([sg_response('CPG001', {'consanguinity': '1'})]) == {'CPG001'}


def test_surrounding_whitespace_is_ignored():
    assert consanguineous_sg_ids([sg_response('CPG001', {'Consanguinity': ' 1 '})]) == {'CPG001'}


def test_anything_other_than_one_fails_closed():
    # A future 'yes' or 'true' must keep the flag rather than silently excusing the pair.
    sgs = [
        sg_response('CPG001', {'Consanguinity': 'yes'}),
        sg_response('CPG002', {'Consanguinity': 'true'}),
        sg_response('CPG003', {'Consanguinity': '2'}),
        sg_response('CPG004', {'Consanguinity': None}),
    ]

    assert consanguineous_sg_ids(sgs) == set()


def test_an_integer_one_is_accepted_too():
    # Recorded as a string in practice, but the phenotypes blob is untyped.
    assert consanguineous_sg_ids([sg_response('CPG001', {'Consanguinity': 1})]) == {'CPG001'}


# ---------------------------------------------------------------------------
# Walking from a flagged co-parent pair to their children
# ---------------------------------------------------------------------------
# CPG001 (dad) and CPG002 (mom) are co-parents of CPG003. Their measured relatedness is
# second-degree, which peddy's 'mom-dad' expectation does not allow.
TRIO_PED = 'FAM1\tCPG001\t0\t0\t1\t1\nFAM1\tCPG002\t0\t0\t2\t1\nFAM1\tCPG003\tCPG001\tCPG002\t1\t1\n'


def related_co_parents(tmp_path: Path, ped: str = TRIO_PED) -> dict[str, str]:
    return write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('CPG001', inferred_sex=MALE, provided_sex='male'),
            sample_row('CPG002', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[],
        ped=ped,
    )


def co_parent_flags(
    inputs: dict[str, str], consanguineous: set[str] | None = None, contaminated: set[str] | None = None
) -> list:
    flags, _, _ = produce_flags(**inputs, consanguineous_sgs=consanguineous, contaminated_sgs=contaminated)
    return [
        flag
        for sg_flags in flags.values()
        for flag in sg_flags
        if getattr(flag, 'expected_relationship', None) == 'mom-dad'
    ]


def with_pair(tmp_path: Path, ped: str = TRIO_PED) -> dict[str, str]:
    """The trio inputs plus the one pairs.tsv row that makes the co-parents look related."""
    from test_check_pedigree import pair_row  # noqa: PLC0415

    return write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('CPG001', inferred_sex=MALE, provided_sex='male'),
            sample_row('CPG002', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[pair_row('CPG001', 'CPG002', relatedness=CONSANGUINEOUS_KIN, ibs0=CONSANGUINEOUS_IBS0)],
        ped=ped,
    )


def test_related_co_parents_are_flagged_when_nothing_is_recorded(tmp_path):
    assert len(co_parent_flags(with_pair(tmp_path))) == 1


def test_a_child_recording_the_union_excuses_the_pair(tmp_path):
    assert co_parent_flags(with_pair(tmp_path), consanguineous={'CPG003'}) == []


def test_an_unrelated_sg_recording_the_union_does_not_excuse_the_pair(tmp_path):
    assert len(co_parent_flags(with_pair(tmp_path), consanguineous={'CPG999'})) == 1


def test_any_child_recording_the_union_is_enough(tmp_path):
    # Real families disagree with themselves: the proband carries the phenotype while a sibling
    # sequenced in the same family records '0'.
    two_kids = TRIO_PED + 'FAM1\tCPG004\tCPG001\tCPG002\t2\t1\n'

    assert co_parent_flags(with_pair(tmp_path, two_kids), consanguineous={'CPG004'}) == []


def test_a_child_of_only_one_parent_does_not_excuse_the_pair(tmp_path):
    # A half sibling says nothing about whether these two particular people are related.
    half_sibling = TRIO_PED + 'FAM1\tCPG005\tCPG001\t0\t1\t1\n'

    assert len(co_parent_flags(with_pair(tmp_path, half_sibling), consanguineous={'CPG005'})) == 1


def test_co_parents_with_no_children_in_the_pedigree_stay_flagged(tmp_path):
    childless = 'FAM1\tCPG001\t0\t0\t1\t1\nFAM1\tCPG002\t0\t0\t2\t1\n'

    # Without a child there is no 'mom-dad' expectation at all, so this is the unrelated case and
    # second-degree relatedness across it remains a conflict.
    flags, _, _ = produce_flags(**with_pair(tmp_path, childless), consanguineous_sgs={'CPG001', 'CPG002'})

    assert [f for sg_flags in flags.values() for f in sg_flags] != []


def test_the_excuse_applies_only_to_co_parents(tmp_path):
    # A parent-child pair measuring unrelated is a swap, and a recorded union cannot explain it.
    from test_check_pedigree import pair_row  # noqa: PLC0415

    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('CPG001', inferred_sex=MALE, provided_sex='male'),
            sample_row('CPG003', inferred_sex=MALE, provided_sex='male'),
        ],
        pair_rows=[pair_row('CPG001', 'CPG003', relatedness=0.01)],
        ped=TRIO_PED,
    )

    flags, _, _ = produce_flags(**inputs, consanguineous_sgs={'CPG003'})

    assert [f for sg_flags in flags.values() for f in sg_flags] != []


def test_omitting_the_argument_keeps_the_previous_behaviour(tmp_path):
    # produce_flags is called without it by anything that has no Metamist access.
    flags, _, _ = produce_flags(**with_pair(tmp_path))

    assert [f for sg_flags in flags.values() for f in sg_flags] != []
