"""
Tests for inferring relatedness from what somalier measured, rather than from its --infer output.

The band edges are calibrated against perth-neuro, where the pedigree already states the
relationship for 99,999 pairs. The `OBSERVED_*` constants below are the real measured ranges from
that dataset, so these tests double as a record of what the calibration was based on: if a band
edge moves far enough to reclassify a real cluster, a test fails.
"""

import pytest

from rd_qc.utils import (
    DEGREE_IDENTICAL,
    DEGREE_PARENT_CHILD,
    DEGREE_SECOND,
    DEGREE_SIBLINGS,
    DEGREE_THIRD,
    DEGREE_UNRELATED,
    FIRST_DEGREE_MIN_RELATEDNESS,
    PARENT_CHILD_MAX_IBS0_RATIO,
    SECOND_DEGREE_MIN_RELATEDNESS,
    THIRD_DEGREE_MIN_RELATEDNESS,
    UNSPECIFIED_RELATED,
    VERDICT_CONFLICT,
    VERDICT_OK,
    VERDICT_REFINEMENT,
    expected_degrees,
    infer_degree,
    refine_expected_relationship,
    relatedness_verdict,
)

# somalier compares ~16.5k sites; the real spread was 14197..17002, only 1.2x.
SITES = 16500

# The measured clusters on perth-neuro, as (kin_min, kin_max, ibs0_min, ibs0_max, n_pairs).
OBSERVED_PARENT_CHILD = (0.438, 0.548, 0, 12, 215)
OBSERVED_FULL_SIBLINGS = (0.430, 0.554, 256, 459, 26)
OBSERVED_GRANDCHILD = (0.203, 0.306, 554, 830, 8)
OBSERVED_NIECE_NEPHEW = (0.204, 0.296, 556, 850, 8)


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
    # Every real pair in these clusters must classify correctly, not just the median.
    kin_min, kin_max, ibs0_min, ibs0_max, _ = cluster

    assert infer_degree(kin_min, ibs0_max, SITES) == expected
    assert infer_degree(kin_max, ibs0_min, SITES) == expected


def test_ibs0_separates_the_two_first_degree_relationships():
    # A parent and child share an allele at every site, so ibs0 is ~0. Full siblings do not.
    # On perth-neuro parent-child reached ibs0 12 and siblings never went below 201.
    assert infer_degree(0.5, 12, SITES) == DEGREE_PARENT_CHILD
    assert infer_degree(0.5, 201, SITES) == DEGREE_SIBLINGS


def test_the_ibs0_threshold_has_margin_on_both_sides():
    observed_parent_child_max = OBSERVED_PARENT_CHILD[3] / SITES
    observed_sibling_min = OBSERVED_FULL_SIBLINGS[2] / SITES

    # The threshold sits strictly between the two observed clusters, with room to spare.
    assert observed_parent_child_max < PARENT_CHILD_MAX_IBS0_RATIO < observed_sibling_min
    assert observed_parent_child_max * 2 < PARENT_CHILD_MAX_IBS0_RATIO
    assert observed_sibling_min / 2 > PARENT_CHILD_MAX_IBS0_RATIO


def test_ibs0_is_normalised_by_the_site_count():
    # The same ibs0 means different things at different site counts, so the threshold is a ratio.
    # 50 out of 16500 sites is 0.3% and reads as parent-child; the same 50 out of 5000 is 1%.
    assert infer_degree(0.5, 50, 16500) == DEGREE_PARENT_CHILD
    assert infer_degree(0.5, 50, 5000) == DEGREE_SIBLINGS


def test_zero_sites_does_not_divide_by_zero():
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
    ],
)
def test_expected_degrees(relationship, acceptable):
    assert expected_degrees(relationship) == acceptable


def test_siblings_accepts_either_half_or_full():
    # peddy says 'siblings' when only one shared parent is recorded, which is satisfied by a half
    # sibling at ~0.25 or a full sibling at ~0.5. This is what kills the biggest false-flag class.
    assert expected_degrees('siblings') == {DEGREE_SIBLINGS, DEGREE_SECOND}
    assert relatedness_verdict('siblings', DEGREE_SIBLINGS) == VERDICT_OK
    assert relatedness_verdict('siblings', DEGREE_SECOND) == VERDICT_OK


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
        # It states one and the measurement contradicts it.
        ('parent-child', DEGREE_UNRELATED, VERDICT_CONFLICT),
        ('parent-child', DEGREE_SIBLINGS, VERDICT_CONFLICT),
        ('unrelated', DEGREE_SIBLINGS, VERDICT_CONFLICT),
        ('mom-dad', DEGREE_SECOND, VERDICT_CONFLICT),
        ('grandchild', DEGREE_PARENT_CHILD, VERDICT_CONFLICT),
        # It records no path between them.
        (UNSPECIFIED_RELATED, DEGREE_UNRELATED, VERDICT_OK),
        (UNSPECIFIED_RELATED, DEGREE_THIRD, VERDICT_REFINEMENT),
        (UNSPECIFIED_RELATED, DEGREE_SECOND, VERDICT_REFINEMENT),
        (UNSPECIFIED_RELATED, DEGREE_SIBLINGS, VERDICT_REFINEMENT),
        (UNSPECIFIED_RELATED, DEGREE_PARENT_CHILD, VERDICT_REFINEMENT),
        (UNSPECIFIED_RELATED, DEGREE_IDENTICAL, VERDICT_CONFLICT),
    ],
)
def test_relatedness_verdict(relationship, measured, verdict):
    assert relatedness_verdict(relationship, measured) == verdict


def test_identical_genomes_are_always_a_conflict():
    # Two samples with one genome is never a pedigree omission: it is one sample recorded twice, or
    # a swap. perth-neuro has exactly one, and the old label-based check missed it entirely.
    for relationship in (UNSPECIFIED_RELATED, 'unknown', 'unrelated', 'parent-child', 'siblings'):
        assert relatedness_verdict(relationship, DEGREE_IDENTICAL) == VERDICT_CONFLICT


def test_in_laws_are_not_flagged():
    # Same family, no recorded path, and the genotypes agree they are unrelated. Nothing to say.
    reframed = refine_expected_relationship('unrelated', 'FAM1', 'FAM1')

    assert relatedness_verdict(reframed, DEGREE_UNRELATED) == VERDICT_OK


def test_the_consanguinity_case_survives():
    # perth-neuro's one real finding: a recorded mother and father measuring second-degree. peddy
    # reports co-parents as 'mom-dad' rather than 'unrelated', so the reframing never touches them.
    assert refine_expected_relationship('mom-dad', 'FAM1', 'FAM1') == 'mom-dad'
    assert infer_degree(0.181, 871, SITES) == DEGREE_SECOND
    assert relatedness_verdict('mom-dad', DEGREE_SECOND) == VERDICT_CONFLICT


# ---------------------------------------------------------------------------
# refine_expected_relationship
# ---------------------------------------------------------------------------
def test_same_family_unrelated_becomes_an_unspecified_expectation():
    assert refine_expected_relationship('unrelated', 'FAM1', 'FAM1') == UNSPECIFIED_RELATED


def test_cross_family_unrelated_is_left_alone():
    # A related measurement across a family boundary must stay a conflict.
    assert refine_expected_relationship('unrelated', 'FAM1', 'FAM2') == 'unrelated'
    assert relatedness_verdict('unrelated', DEGREE_SECOND) == VERDICT_CONFLICT


def test_unrelated_with_an_unknown_family_is_left_alone():
    assert refine_expected_relationship('unrelated', None, None) == 'unrelated'
    assert refine_expected_relationship('unrelated', '', 'FAM1') == 'unrelated'


@pytest.mark.parametrize(
    'relationship',
    ['parent-child', 'full siblings', 'siblings', 'grandchild', 'niece/nephew', 'cousins', 'mom-dad'],
)
def test_every_other_relationship_passes_through_untouched(relationship):
    assert refine_expected_relationship(relationship, 'FAM1', 'FAM1') == relationship


# ---------------------------------------------------------------------------
# Calibration against the measured background
# ---------------------------------------------------------------------------
# Cross-family pairs on perth-neuro, which the pedigree really does expect to be unrelated, over
# 99,515 pairs. No upper cluster, just a smooth tail, so anything inside it is background.
OBSERVED_CROSS_FAMILY_MEDIAN = -0.006
OBSERVED_CROSS_FAMILY_P99_9 = 0.072
OBSERVED_CROSS_FAMILY_MAX = 0.158
# The nearest real second-degree cluster starts here, leaving a genuine gap above the background.
OBSERVED_SECOND_DEGREE_MIN = 0.203


def test_the_second_degree_bound_sits_in_the_gap_above_background():
    assert OBSERVED_CROSS_FAMILY_MAX < SECOND_DEGREE_MIN_RELATEDNESS < OBSERVED_SECOND_DEGREE_MIN


def test_the_measured_background_never_reaches_second_degree():
    # Every one of those 99,515 pairs must read as third-degree or unrelated, never closer,
    # otherwise the report claims a pedigree error where there is only cohort background.
    for kin in (OBSERVED_CROSS_FAMILY_MEDIAN, OBSERVED_CROSS_FAMILY_P99_9, OBSERVED_CROSS_FAMILY_MAX):
        assert infer_degree(kin, 1400, SITES) in (DEGREE_UNRELATED, DEGREE_THIRD)


def test_bounds_are_the_geometric_midpoints_between_degrees():
    # somalier reports 2*phi, so successive degrees are 1.0, 0.5, 0.25, 0.125 and the split
    # between two adjacent expectations is their geometric mean.
    assert pytest.approx((0.5 * 0.25) ** 0.5, abs=0.005) == FIRST_DEGREE_MIN_RELATEDNESS
    assert pytest.approx((0.25 * 0.125) ** 0.5, abs=0.005) == SECOND_DEGREE_MIN_RELATEDNESS
    assert pytest.approx((0.125 * 0.0625) ** 0.5, abs=0.005) == THIRD_DEGREE_MIN_RELATEDNESS


def test_distant_relatedness_against_an_unrelated_expectation_is_only_a_refinement():
    # Not assertable: a real first cousin sits at 0.125, and this cohort's background p99.99 was
    # 0.126. The 40 pairs it caught were concentrated on a few samples, one in 8 pairs each.
    assert relatedness_verdict('unrelated', DEGREE_THIRD) == VERDICT_REFINEMENT
    assert relatedness_verdict('mom-dad', DEGREE_THIRD) == VERDICT_REFINEMENT


def test_closer_than_third_degree_against_an_unrelated_expectation_is_still_a_conflict():
    # Second-degree and above are clear of the background, so they remain assertable.
    assert relatedness_verdict('unrelated', DEGREE_SECOND) == VERDICT_CONFLICT
    assert relatedness_verdict('unrelated', DEGREE_SIBLINGS) == VERDICT_CONFLICT
    assert relatedness_verdict('unrelated', DEGREE_IDENTICAL) == VERDICT_CONFLICT


def test_a_cousin_expectation_is_still_satisfied_by_a_third_degree_measurement():
    # The demotion above is specific to an 'unrelated' expectation; a stated cousin relationship
    # confirmed by the measurement must not become a refinement.
    assert relatedness_verdict('cousins', DEGREE_THIRD) == VERDICT_OK
