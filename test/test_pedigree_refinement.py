"""
Tests for the conflict/refinement classification shared by the pedigree check and the report.

The cases below are the complete expected->inferred transition table observed on perth-neuro
(163 pedigree flags across 480 SGs), so this doubles as a regression record of what real data
actually produces. 92 of those 163 are refinements.
"""

import pytest
from peddy import Ped

from rd_qc.utils import (
    PEDDY_RELATIONSHIPS,
    UNSPECIFIED_RELATED,
    is_pedigree_refinement,
    refine_expected_relationship,
)

# (expected, inferred, count observed on perth-neuro)
REFINEMENTS = [
    ('siblings', 'full siblings', 37),
    (UNSPECIFIED_RELATED, 'niece/nephew', 36),
    (UNSPECIFIED_RELATED, 'parent-child', 10),
    (UNSPECIFIED_RELATED, 'full siblings', 5),
    (UNSPECIFIED_RELATED, 'grandchild', 4),
]

CONFLICTS = [
    ('unrelated', UNSPECIFIED_RELATED, 26),
    ('unrelated', 'full siblings', 21),
    ('unrelated', 'niece/nephew', 10),
    ('unrelated', 'parent-child', 7),
    ('parent-child', 'unrelated', 2),
    ('unrelated', 'mom-dad', 2),
    ('grandchild', UNSPECIFIED_RELATED, 1),
    ('mom-dad', UNSPECIFIED_RELATED, 1),
    ('parent-child', UNSPECIFIED_RELATED, 1),
]

REFINEMENT_TOTAL = 92
CONFLICT_TOTAL = 71


@pytest.mark.parametrize(('expected', 'inferred'), [(e, i) for e, i, _ in REFINEMENTS])
def test_refinements_are_classified_as_refinements(expected, inferred):
    assert is_pedigree_refinement(expected, inferred) is True


@pytest.mark.parametrize(('expected', 'inferred'), [(e, i) for e, i, _ in CONFLICTS])
def test_conflicts_are_not_classified_as_refinements(expected, inferred):
    assert is_pedigree_refinement(expected, inferred) is False


def test_the_observed_transition_table_splits_92_to_71():
    assert sum(n for _, _, n in REFINEMENTS) == REFINEMENT_TOTAL
    assert sum(n for _, _, n in CONFLICTS) == CONFLICT_TOTAL


def test_an_unspecified_expected_relationship_contradicted_by_unrelated_is_a_conflict():
    # The pedigree says these two are related somehow; the genotypes say they are not. That is a
    # contradiction, not the pedigree being vague.
    assert is_pedigree_refinement(UNSPECIFIED_RELATED, 'unrelated') is False


def test_full_siblings_downgraded_to_siblings_is_a_conflict():
    # The reverse of the common case: the pedigree claimed both parents are shared and the
    # genotypes only support one. Worth a look.
    assert is_pedigree_refinement('full siblings', 'siblings') is False


def test_a_missing_pedigree_entry_is_a_conflict():
    # check_pedigree falls back to 'unknown' when a sample is absent from the PED entirely. That
    # is missing data someone should fix, so it stays visible rather than being demoted.
    assert is_pedigree_refinement('unknown', 'parent-child') is False


def test_peddy_vocabulary_is_pinned():
    # The classification compares these strings literally, so an upstream rename must break a test
    # rather than silently reclassifying every flag in the dataset.
    assert UNSPECIFIED_RELATED in PEDDY_RELATIONSHIPS
    assert {'siblings', 'full siblings', 'unrelated', 'parent-child'} <= PEDDY_RELATIONSHIPS
    assert hasattr(Ped, 'relation')


@pytest.mark.parametrize('expected', sorted(PEDDY_RELATIONSHIPS))
@pytest.mark.parametrize('inferred', sorted(PEDDY_RELATIONSHIPS))
def test_classification_is_total_over_the_vocabulary(expected, inferred):
    # Never raises, and always returns a bool, for every pair peddy can produce.
    assert isinstance(is_pedigree_refinement(expected, inferred), bool)


# ---------------------------------------------------------------------------
# Reinterpreting peddy's 'unrelated' within a family
# ---------------------------------------------------------------------------
def test_same_family_unrelated_becomes_related_at_unknown_level():
    # peddy says 'unrelated' for two family members with no recorded blood path. All 66 of
    # perth-neuro's 'expected unrelated' flags were this case, and none were cross-family.
    assert refine_expected_relationship('unrelated', 'full siblings', 'FAM1', 'FAM1') == UNSPECIFIED_RELATED


def test_cross_family_unrelated_is_left_alone():
    # Genuinely unrelated, and a related inference here is the cross-family swap case.
    assert refine_expected_relationship('unrelated', 'full siblings', 'FAM1', 'FAM2') == 'unrelated'


def test_unrelated_with_an_unknown_family_is_left_alone():
    assert refine_expected_relationship('unrelated', 'full siblings', None, None) == 'unrelated'
    assert refine_expected_relationship('unrelated', 'full siblings', '', 'FAM1') == 'unrelated'


@pytest.mark.parametrize('relation', sorted(PEDDY_RELATIONSHIPS - {'unrelated'}))
def test_every_other_relationship_passes_through_untouched(relation):
    # Only 'unrelated' is reinterpreted. In particular 'mom-dad' survives, which is what keeps
    # consanguinity between two recorded parents visible as a conflict.
    assert refine_expected_relationship(relation, 'full siblings', 'FAM1', 'FAM1') == relation


def test_same_family_reframing_turns_a_false_alarm_into_a_refinement():
    # The end-to-end effect: 'unrelated' -> 'full siblings' within one family was a conflict, and
    # is now correctly a refinement, because the pedigree never claimed they were unrelated.
    reframed = refine_expected_relationship('unrelated', 'full siblings', 'FAM1', 'FAM1')

    assert is_pedigree_refinement('unrelated', 'full siblings') is False
    assert is_pedigree_refinement(reframed, 'full siblings') is True


def test_reframing_never_creates_a_mismatch_out_of_an_agreeing_pair():
    """
    Two family members with no recorded link whom somalier also calls unrelated agree, and must
    keep agreeing. Reframing them to 'related at unknown level' would invent a conflict: on
    perth-neuro that was 18 new conflicts drawn from pairs that previously matched, which is the
    in-law case arriving as an alarm.
    """
    assert refine_expected_relationship('unrelated', 'unrelated', 'FAM1', 'FAM1') == 'unrelated'
