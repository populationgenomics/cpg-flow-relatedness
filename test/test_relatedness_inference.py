"""
Tests for inferring relatedness from what somalier measured, rather than from its --infer output.

The `OBSERVED_*` constants are measured kinship ranges from a calibration cohort of ~100k pairs
whose relationships the pedigree already states. They pin the band edges to real data: if an edge
moves far enough to reclassify a known cluster, a test fails. An edge moving *within* the gaps
between clusters reclassifies nothing a caller can observe, and deliberately fails nothing.
"""

from typing import NamedTuple

import pytest

from rd_qc.utils import (
    DEGREE_IDENTICAL,
    DEGREE_PARENT_CHILD,
    DEGREE_SECOND,
    DEGREE_SIBLINGS,
    DEGREE_THIRD,
    DEGREE_UNRELATED,
    PARENT_CHILD_MAX_IBS0_RATIO,
    SECOND_DEGREE_MIN_RELATEDNESS,
    UNSPECIFIED_RELATED,
    VERDICT_CONFLICT,
    VERDICT_OK,
    VERDICT_REFINEMENT,
    expected_degrees,
    infer_degree,
    refine_expected_relationship,
    relatedness_verdict,
)

# somalier compares ~16.5k sites, with little spread between pairs.
SITES = 16500


class Cluster(NamedTuple):
    """One measured relationship's observed range, so a test can name the edge it compares."""

    kin_min: float
    kin_max: float
    ibs0_min: int
    ibs0_max: int
    n_pairs: int


OBSERVED_PARENT_CHILD = Cluster(0.438, 0.548, 0, 12, 215)
OBSERVED_FULL_SIBLINGS = Cluster(0.430, 0.554, 256, 459, 26)
OBSERVED_GRANDCHILD = Cluster(0.203, 0.306, 554, 830, 8)
OBSERVED_NIECE_NEPHEW = Cluster(0.204, 0.296, 556, 850, 8)


# ---------------------------------------------------------------------------
# infer_degree
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ('relatedness', 'ibs0', 'expected'),
    [
        (1.000, 0, DEGREE_IDENTICAL),
        (0.95, 0, DEGREE_IDENTICAL),
        (0.50, 0, DEGREE_PARENT_CHILD),
        (0.50, 340, DEGREE_SIBLINGS),
        (0.25, 700, DEGREE_SECOND),
        (0.12, 1100, DEGREE_THIRD),
        (0.00, 1414, DEGREE_UNRELATED),
        (-0.577, 1400, DEGREE_UNRELATED),
    ],
)
def test_infer_degree(relatedness, ibs0, expected):
    assert infer_degree(relatedness, ibs0, SITES) == expected


@pytest.mark.parametrize(
    ('cluster', 'expected'),
    [
        (OBSERVED_PARENT_CHILD, DEGREE_PARENT_CHILD),
        (OBSERVED_FULL_SIBLINGS, DEGREE_SIBLINGS),
        (OBSERVED_GRANDCHILD, DEGREE_SECOND),
        (OBSERVED_NIECE_NEPHEW, DEGREE_SECOND),
    ],
)
def test_observed_clusters_land_in_the_right_band_at_both_extremes(cluster, expected):
    # Both extremes of each cluster must classify correctly, not just the median.
    assert infer_degree(cluster.kin_min, cluster.ibs0_max, SITES) == expected
    assert infer_degree(cluster.kin_max, cluster.ibs0_min, SITES) == expected


def test_ibs0_separates_the_two_first_degree_relationships():
    # A parent and child share an allele at every site, so ibs0 is ~0. Full siblings do not.
    assert infer_degree(0.5, 12, SITES) == DEGREE_PARENT_CHILD
    assert infer_degree(0.5, 201, SITES) == DEGREE_SIBLINGS


def test_the_ibs0_threshold_sits_between_the_two_observed_clusters():
    # The only property that matters: no observed parent-child pair reads as siblings and no
    # observed sibling pair reads as parent-child. Where inside the gap the edge sits is a tuning
    # decision no caller can observe.
    assert (
        OBSERVED_PARENT_CHILD.ibs0_max / SITES < PARENT_CHILD_MAX_IBS0_RATIO < OBSERVED_FULL_SIBLINGS.ibs0_min / SITES
    )


def test_ibs0_is_normalised_by_the_site_count():
    # The same ibs0 means different things at different site counts, so the threshold is a ratio.
    # 50 of 16500 sites is 0.3% and reads as parent-child; the same 50 of 500 is 10%.
    assert infer_degree(0.5, 50, 16500) == DEGREE_PARENT_CHILD
    assert infer_degree(0.5, 50, 500) == DEGREE_SIBLINGS


def test_a_pair_with_no_compared_sites_reads_as_parent_child():
    # No sites compared means no ibs0 evidence either way, so the ratio floors to 0.0 and the
    # kinship alone decides. Documented because the alternative was a ZeroDivisionError.
    assert infer_degree(0.5, 0, 0) == DEGREE_PARENT_CHILD


# ---------------------------------------------------------------------------
# expected_degrees
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ('relationship', 'acceptable'),
    [
        ('parent-child', {DEGREE_PARENT_CHILD}),
        ('full siblings', {DEGREE_SIBLINGS}),
        ('grandchild', {DEGREE_SECOND}),
        ('niece/nephew', {DEGREE_SECOND}),
        ('cousins', {DEGREE_THIRD}),
        ('mom-dad', {DEGREE_UNRELATED}),
        ('unrelated', {DEGREE_UNRELATED}),
        # peddy says 'siblings' when only one shared parent is recorded, which is satisfied by a
        # half sibling at ~0.25 or a full sibling at ~0.5.
        ('siblings', {DEGREE_SIBLINGS, DEGREE_SECOND}),
    ],
)
def test_expected_degrees(relationship, acceptable):
    assert expected_degrees(relationship) == acceptable


@pytest.mark.parametrize('relationship', [UNSPECIFIED_RELATED, 'unknown', 'something new'])
def test_an_unusable_expectation_has_no_acceptable_degrees(relationship):
    assert expected_degrees(relationship) is None


# ---------------------------------------------------------------------------
# relatedness_verdict
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ('relationship', 'measured', 'verdict'),
    [
        # The pedigree states a relationship and the measurement agrees.
        ('parent-child', DEGREE_PARENT_CHILD, VERDICT_OK),
        ('full siblings', DEGREE_SIBLINGS, VERDICT_OK),
        ('grandchild', DEGREE_SECOND, VERDICT_OK),
        ('unrelated', DEGREE_UNRELATED, VERDICT_OK),
        ('mom-dad', DEGREE_UNRELATED, VERDICT_OK),
        # A half-sibling expectation is satisfied at either degree.
        ('siblings', DEGREE_SIBLINGS, VERDICT_OK),
        ('siblings', DEGREE_SECOND, VERDICT_OK),
        # A stated cousin expectation confirmed by the measurement. Reads OK for a different
        # reason than the cross-family third-degree row below: here the pedigree said so.
        ('cousins', DEGREE_THIRD, VERDICT_OK),
        # It states one and the measurement contradicts it.
        ('parent-child', DEGREE_UNRELATED, VERDICT_CONFLICT),
        ('parent-child', DEGREE_SIBLINGS, VERDICT_CONFLICT),
        ('unrelated', DEGREE_SIBLINGS, VERDICT_CONFLICT),
        ('unrelated', DEGREE_SECOND, VERDICT_CONFLICT),
        ('mom-dad', DEGREE_SECOND, VERDICT_CONFLICT),
        ('grandchild', DEGREE_PARENT_CHILD, VERDICT_CONFLICT),
        # Third-degree across a family boundary is not assertable: a real first cousin sits at
        # 0.125, inside the cohort background tail, so it is indistinguishable from background.
        # An 'unrelated' expectation only survives refinement when the pair is cross-family.
        ('unrelated', DEGREE_THIRD, VERDICT_OK),
        # Same measurement, but co-parents are within one family, where distant relatedness
        # speaks to consanguinity in that family rather than to cohort background.
        ('mom-dad', DEGREE_THIRD, VERDICT_REFINEMENT),
        # It records no path between them.
        (UNSPECIFIED_RELATED, DEGREE_UNRELATED, VERDICT_OK),
        (UNSPECIFIED_RELATED, DEGREE_THIRD, VERDICT_REFINEMENT),
        (UNSPECIFIED_RELATED, DEGREE_SECOND, VERDICT_REFINEMENT),
        (UNSPECIFIED_RELATED, DEGREE_SIBLINGS, VERDICT_REFINEMENT),
        (UNSPECIFIED_RELATED, DEGREE_PARENT_CHILD, VERDICT_REFINEMENT),
        # Two samples with one genome is never a pedigree omission: it is one sample recorded
        # twice, or a swap. So it is a conflict against every expectation, stated or not.
        (UNSPECIFIED_RELATED, DEGREE_IDENTICAL, VERDICT_CONFLICT),
        ('unknown', DEGREE_IDENTICAL, VERDICT_CONFLICT),
        ('unrelated', DEGREE_IDENTICAL, VERDICT_CONFLICT),
        ('parent-child', DEGREE_IDENTICAL, VERDICT_CONFLICT),
        ('siblings', DEGREE_IDENTICAL, VERDICT_CONFLICT),
    ],
)
def test_relatedness_verdict(relationship, measured, verdict):
    assert relatedness_verdict(relationship, measured) == verdict


# ---------------------------------------------------------------------------
# refine_expected_relationship
# ---------------------------------------------------------------------------
def test_same_family_unrelated_becomes_an_unspecified_expectation():
    # Two family members with no recorded blood path between them. peddy calls that 'unrelated',
    # but the pedigree never asserted it, so the expectation is reframed rather than trusted.
    # This is the coupling that makes the cross-family third-degree silence above safe: a
    # same-family pair never reaches relatedness_verdict still carrying 'unrelated'.
    assert refine_expected_relationship('unrelated', 'FAM1', 'FAM1') == UNSPECIFIED_RELATED


def test_cross_family_unrelated_is_left_alone():
    # Two people in different families really are expected to be unrelated, so a related
    # measurement across the boundary has to stay assertable.
    assert refine_expected_relationship('unrelated', 'FAM1', 'FAM2') == 'unrelated'


def test_unrelated_with_an_unknown_family_is_left_alone():
    assert refine_expected_relationship('unrelated', None, None) == 'unrelated'
    assert refine_expected_relationship('unrelated', '', 'FAM1') == 'unrelated'


@pytest.mark.parametrize(
    'relationship',
    ['parent-child', 'full siblings', 'siblings', 'grandchild', 'niece/nephew', 'cousins', 'mom-dad'],
)
def test_every_other_relationship_passes_through_untouched(relationship):
    # peddy reports co-parents as 'mom-dad' rather than 'unrelated', so the reframing never
    # touches them and a consanguineous union stays surfaced.
    assert refine_expected_relationship(relationship, 'FAM1', 'FAM1') == relationship


# ---------------------------------------------------------------------------
# Calibration against the measured background
# ---------------------------------------------------------------------------
# Cross-family pairs, which the pedigree really does expect to be unrelated. These form no upper
# cluster, just a smooth tail, so anything inside the tail is cohort background.
OBSERVED_CROSS_FAMILY_MEDIAN = -0.006
OBSERVED_CROSS_FAMILY_P99_9 = 0.072
OBSERVED_CROSS_FAMILY_MAX = 0.158
# The nearest real second-degree cluster starts here, leaving a genuine gap above the background.
OBSERVED_SECOND_DEGREE_MIN = 0.203


def test_the_second_degree_bound_sits_in_the_gap_above_background():
    assert OBSERVED_CROSS_FAMILY_MAX < SECOND_DEGREE_MIN_RELATEDNESS < OBSERVED_SECOND_DEGREE_MIN


@pytest.mark.parametrize(
    ('kin', 'expected'),
    [
        (OBSERVED_CROSS_FAMILY_MEDIAN, DEGREE_UNRELATED),
        (OBSERVED_CROSS_FAMILY_P99_9, DEGREE_UNRELATED),
        (OBSERVED_CROSS_FAMILY_MAX, DEGREE_THIRD),
    ],
)
def test_the_measured_background_never_reaches_second_degree(kin, expected):
    # The whole background tail must read as third-degree or unrelated, never closer, otherwise
    # the report claims a pedigree error where there is only cohort background.
    assert infer_degree(kin, 1400, SITES) == expected
