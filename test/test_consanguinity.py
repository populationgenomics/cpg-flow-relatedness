"""
Tests for excusing a co-parent relatedness flag with a recorded consanguineous union.

Metamist records consanguinity against the child, in the free-text participant phenotypes, so the
check has to walk from a flagged mom-dad pair to their children. Two halves are tested here: the
extraction of the phenotype out of a Metamist response, and the pedigree walk inside produce_flags.

Absence is deliberately treated as "not recorded", so a dataset that has never supplied the field
keeps every co-parent flag rather than having them all quietly excused.
"""

from pathlib import Path

import pytest
from fixtures.somalier_tsv import (
    FEMALE,
    MALE,
    TRIO_PED,
    pair_row,
    sample_row,
    write_inputs,
)

from rd_qc.scripts.check_pedigree import CO_PARENTS, produce_flags
from rd_qc.utils import (
    DEGREE_SECOND,
    DEGREE_UNRELATED,
    UNSPECIFIED_RELATED,
    VERDICT_REFINEMENT,
    SomalierRelatednessFlag,
    consanguineous_sg_ids,
)

# Second-degree relatedness: above SECOND_DEGREE_MIN_RELATEDNESS, below first-degree.
CONSANGUINEOUS_KIN = 0.21
CONSANGUINEOUS_IBS0 = 800


def sg_response(sg_id: str, phenotypes: dict | None) -> dict:
    """One sequencingGroups entry in the shape SG_QUERY returns."""
    return {'id': sg_id, 'sample': {'participant': {'phenotypes': phenotypes}}}


# ---------------------------------------------------------------------------
# Reading the phenotype out of Metamist's response
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ('phenotypes', 'recorded'),
    [
        # Only a literal 1 records a union. Anything else fails closed, keeping the flag.
        ({'Consanguinity': '1'}, {'TST001'}),
        # Recorded as a string in practice, but the phenotypes blob is untyped.
        ({'Consanguinity': 1}, {'TST001'}),
        # 'Consanguinity' is the house spelling, but the phenotypes dict is free text.
        ({'consanguinity': '1'}, {'TST001'}),
        ({'Consanguinity': ' 1 '}, {'TST001'}),
        ({'Consanguinity': '0'}, set()),
        # A future 'yes' or 'true' must keep the flag rather than silently excusing the pair.
        ({'Consanguinity': 'yes'}, set()),
        ({'Consanguinity': 'true'}, set()),
        ({'Consanguinity': '2'}, set()),
        ({'Consanguinity': None}, set()),
        ({'Birth Year': '1984'}, set()),
        ({}, set()),
        (None, set()),
    ],
)
def test_only_a_literal_one_records_a_consanguineous_union(phenotypes, recorded):
    assert consanguineous_sg_ids([sg_response('TST001', phenotypes)]) == recorded


def test_an_sg_with_no_sample_or_participant_records_nothing():
    # Metamist omits the whole branch rather than nesting empty dicts, so neither level is safe
    # to assume. Raising here would take down the pedigree check for the whole dataset.
    assert consanguineous_sg_ids([{'id': 'TST001'}, {'id': 'TST002', 'sample': None}]) == set()


# ---------------------------------------------------------------------------
# Walking from a flagged co-parent pair to their children
# ---------------------------------------------------------------------------
def co_parent_flags(inputs: dict[str, str], consanguineous: set[str] | None = None) -> list:
    """Just the co-parent relatedness flags `produce_flags` raised for these inputs."""
    flags, _, _ = produce_flags(**inputs, consanguineous_sgs=consanguineous)
    return [
        flag
        for sg_flags in flags.values()
        for flag in sg_flags
        if isinstance(flag, SomalierRelatednessFlag) and flag.expected_relationship == CO_PARENTS
    ]


def with_pair(tmp_path: Path, ped: str = TRIO_PED) -> dict[str, str]:
    """The trio inputs plus the one pairs.tsv row that makes the co-parents look related."""
    return write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST002', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[pair_row('TST001', 'TST002', relatedness=CONSANGUINEOUS_KIN, ibs0=CONSANGUINEOUS_IBS0)],
        ped=ped,
    )


def test_related_co_parents_are_flagged_when_nothing_is_recorded(tmp_path):
    assert len(co_parent_flags(with_pair(tmp_path))) == 1


def test_a_child_recording_the_union_excuses_the_pair(tmp_path):
    assert co_parent_flags(with_pair(tmp_path), consanguineous={'TST003'}) == []


def test_an_unrelated_sg_recording_the_union_does_not_excuse_the_pair(tmp_path):
    assert len(co_parent_flags(with_pair(tmp_path), consanguineous={'TST999'})) == 1


def test_any_child_recording_the_union_is_enough(tmp_path):
    # Real families disagree with themselves: the proband carries the phenotype while a sibling
    # sequenced in the same family records '0'.
    two_kids = TRIO_PED + 'FAM1\tTST004\tTST001\tTST002\t2\t1\n'

    assert co_parent_flags(with_pair(tmp_path, two_kids), consanguineous={'TST004'}) == []


def test_a_child_of_only_one_parent_does_not_excuse_the_pair(tmp_path):
    # A half sibling says nothing about whether these two particular people are related.
    half_sibling = TRIO_PED + 'FAM1\tTST005\tTST001\t0\t1\t1\n'

    assert len(co_parent_flags(with_pair(tmp_path, half_sibling), consanguineous={'TST005'})) == 1


def test_co_parents_with_no_children_land_as_a_refinement(tmp_path):
    # Without a child there is no 'mom-dad' expectation at all: peddy calls them 'unrelated',
    # which the same-family reframing turns into "no relationship recorded". So the pair is
    # surfaced as a refinement rather than a conflict, and the recorded union is never consulted.
    # Only-parents-sequenced is rare enough that the de-emphasised section is the right home.
    childless = 'FAM1\tTST001\t0\t0\t1\t1\nFAM1\tTST002\t0\t0\t2\t1\n'

    flags, _, _ = produce_flags(**with_pair(tmp_path, childless), consanguineous_sgs={'TST001', 'TST002'})

    flag = next(f for fs in flags.values() for f in fs)
    assert (flag.expected_relationship, flag.inferred_relationship, flag.verdict) == (
        UNSPECIFIED_RELATED,
        DEGREE_SECOND,
        VERDICT_REFINEMENT,
    )


def test_the_excuse_applies_only_to_co_parents(tmp_path):
    # A parent-child pair measuring unrelated is a swap, and a recorded union cannot explain it.
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('TST001', inferred_sex=MALE, provided_sex='male'),
            sample_row('TST003', inferred_sex=MALE, provided_sex='male'),
        ],
        pair_rows=[pair_row('TST001', 'TST003', relatedness=0.01)],
        ped=TRIO_PED,
    )

    flags, _, _ = produce_flags(**inputs, consanguineous_sgs={'TST003'})

    flagged = [(f.expected_relationship, f.inferred_relationship) for fs in flags.values() for f in fs]

    assert flagged == [('parent-child', DEGREE_UNRELATED)]
