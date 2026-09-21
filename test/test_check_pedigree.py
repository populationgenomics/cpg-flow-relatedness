"""
Tests for the flag-producing half of the pedigree check.

`produce_flags` is the part that reads somalier's relate outputs and decides what is wrong. It
touches neither Metamist nor Slack, which is what lets the local report harness reuse it, so these
tests drive it directly with synthetic somalier TSVs written to a tmp_path.

The column lists below mirror a real somalier 0.3.1 samples.tsv / pairs.tsv header exactly, because
`_check_sex` indexes into them by name and peddy re-reads samples.tsv as a PED file.
"""

from pathlib import Path

from rd_qc.scripts.check_pedigree import produce_flags

SAMPLE_COLUMNS = [
    '#family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype',
    'original_pedigree_sex', 'gt_depth_mean', 'gt_depth_sd', 'depth_mean', 'depth_sd',
    'ab_mean', 'ab_std', 'n_hom_ref', 'n_het', 'n_hom_alt', 'n_unknown', 'p_middling_ab',
    'X_depth_mean', 'X_n', 'X_hom_ref', 'X_het', 'X_hom_alt', 'Y_depth_mean', 'Y_n',
]  # fmt: skip

PAIR_COLUMNS = [
    '#sample_a', 'sample_b', 'relatedness', 'ibs0', 'ibs2', 'hom_concordance', 'hets_a',
    'hets_b', 'hets_ab', 'shared_hets', 'hom_alts_a', 'hom_alts_b', 'shared_hom_alts', 'n',
    'x_ibs0', 'x_ibs2', 'expected_relatedness',
]  # fmt: skip

MALE = 1
FEMALE = 2


def sample_row(
    sample_id: str,
    inferred_sex: int,
    provided_sex: str,
    family: str = 'FAM1',
    paternal: str = '0',
    maternal: str = '0',
    gt_depth_mean: float = 30.0,
) -> str:
    """One somalier samples.tsv row. `inferred_sex` is somalier's call, 1 male / 2 female."""
    values = {
        '#family_id': family, 'sample_id': sample_id, 'paternal_id': paternal,
        'maternal_id': maternal, 'sex': inferred_sex, 'phenotype': 1,
        'original_pedigree_sex': provided_sex, 'gt_depth_mean': gt_depth_mean,
        'gt_depth_sd': 2.0, 'depth_mean': 30.0, 'depth_sd': 2.0, 'ab_mean': 0.5, 'ab_std': 0.1,
        'n_hom_ref': 100, 'n_het': 50, 'n_hom_alt': 40, 'n_unknown': 0, 'p_middling_ab': 0.01,
        'X_depth_mean': 15.0, 'X_n': 500, 'X_hom_ref': 200, 'X_het': 5, 'X_hom_alt': 100,
        'Y_depth_mean': 15.0, 'Y_n': 200,
    }  # fmt: skip
    return '\t'.join(str(values[column]) for column in SAMPLE_COLUMNS)


# somalier compares ~16.5k sites, and ibs0 is judged as a fraction of that, so tests must use a
# realistic count for the parent-child / siblings split to behave as it does in production.
SITES = 16500
IBS0_PARENT_CHILD = 0
IBS0_SIBLINGS = 340
IBS0_UNRELATED = 1414


def pair_row(sample_a: str, sample_b: str, relatedness: float, ibs0: int = IBS0_UNRELATED) -> str:
    """One somalier pairs.tsv row. relatedness and ibs0 are what the check now infers from."""
    values = {
        '#sample_a': sample_a, 'sample_b': sample_b, 'relatedness': relatedness, 'ibs0': ibs0,
        'ibs2': 2000, 'hom_concordance': 0.9, 'hets_a': 50, 'hets_b': 50, 'hets_ab': 60,
        'shared_hets': 40, 'hom_alts_a': 40, 'hom_alts_b': 40, 'shared_hom_alts': 30, 'n': SITES,
        'x_ibs0': 1, 'x_ibs2': 10, 'expected_relatedness': -1,
    }  # fmt: skip
    return '\t'.join(str(values[column]) for column in PAIR_COLUMNS)


def write_inputs(tmp_path: Path, sample_rows: list[str], pair_rows: list[str], ped: str) -> dict[str, str]:
    samples = tmp_path / 'x.samples.tsv'
    pairs = tmp_path / 'x.pairs.tsv'
    expected_ped = tmp_path / 'x.expected.ped'
    samples.write_text('\t'.join(SAMPLE_COLUMNS) + '\n' + '\n'.join(sample_rows) + '\n')
    pairs.write_text('\t'.join(PAIR_COLUMNS) + '\n' + '\n'.join(pair_rows) + '\n')
    expected_ped.write_text(ped)
    return {
        'somalier_samples': str(samples),
        'somalier_pairs': str(pairs),
        'expected_ped_path': str(expected_ped),
    }


# A trio where the expected PED says CPG003 is the child of CPG001 and CPG002, but the measured
# relatedness says all three are unrelated. The two parent-child pairs therefore conflict. The
# CPG001/CPG002 pair does not: peddy calls co-parents 'mom-dad', which expects unrelated anyway.
BROKEN_TRIO_PED = 'FAM1\tCPG001\t0\t0\t1\t1\nFAM1\tCPG002\t0\t0\t2\t1\nFAM1\tCPG003\tCPG001\tCPG002\t1\t1\n'


def broken_trio(tmp_path: Path, provided_sex_003: str = 'male') -> dict[str, str]:
    return write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('CPG001', inferred_sex=MALE, provided_sex='male'),
            sample_row('CPG002', inferred_sex=FEMALE, provided_sex='female'),
            sample_row('CPG003', inferred_sex=MALE, provided_sex=provided_sex_003),
        ],
        pair_rows=[
            pair_row('CPG001', 'CPG002', relatedness=0.01),
            pair_row('CPG001', 'CPG003', relatedness=0.02),
            pair_row('CPG002', 'CPG003', relatedness=0.03),
        ],
        ped=BROKEN_TRIO_PED,
    )


def test_matching_pedigree_and_sex_produces_no_flags(tmp_path):
    # The measurements confirm the PED: both parent-child pairs at ~0.5 with ibs0 ~0, the two
    # co-parents unrelated, and every sex agreeing.
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('CPG001', inferred_sex=MALE, provided_sex='male'),
            sample_row('CPG002', inferred_sex=FEMALE, provided_sex='female'),
            sample_row('CPG003', inferred_sex=MALE, provided_sex='male'),
        ],
        pair_rows=[
            pair_row('CPG001', 'CPG002', relatedness=0.01),
            pair_row('CPG001', 'CPG003', relatedness=0.49, ibs0=IBS0_PARENT_CHILD),
            pair_row('CPG002', 'CPG003', relatedness=0.51, ibs0=IBS0_PARENT_CHILD),
        ],
        ped=BROKEN_TRIO_PED,
    )

    flags, _, _ = produce_flags(**inputs)

    assert flags == {}


def test_sex_mismatch_is_flagged_against_the_owning_sg(tmp_path):
    flags, _, _ = produce_flags(**broken_trio(tmp_path, provided_sex_003='female'))

    sex_flags = [f for f in flags.get('CPG003', []) if f.category == 'sex_inference_mismatch']

    assert len(sex_flags) == 1
    assert (sex_flags[0].provided, sex_flags[0].inferred) == ('female', 'male')


def test_matching_sex_is_not_flagged(tmp_path):
    flags, _, _ = produce_flags(**broken_trio(tmp_path, provided_sex_003='male'))

    assert not [f for fs in flags.values() for f in fs if f.category == 'sex_inference_mismatch']


def test_pedigree_mismatch_is_flagged_for_the_disagreeing_pairs(tmp_path):
    flags, _, _ = produce_flags(**broken_trio(tmp_path))

    pedigree = [f for fs in flags.values() for f in fs if f.category == 'relatedness_mismatch']
    pairs = {(f.sg_id_1, f.sg_id_2) for f in pedigree}

    # Both parent-child pairs conflict. CPG001/CPG002 are the co-parents, whom peddy calls
    # 'mom-dad' and which expects unrelated, so the measurement agrees and raises nothing.
    assert pairs == {('CPG001', 'CPG003'), ('CPG002', 'CPG003')}
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
        ('CPG001', 'CPG003'): 'CPG001',
        ('CPG002', 'CPG003'): 'CPG002',
    }
    assert 'CPG003' not in flags


def test_zero_depth_samples_are_excluded_from_both_checks(tmp_path):
    # A sample somalier could not genotype: gt_depth_mean of 0 means its calls are meaningless.
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('CPG001', inferred_sex=MALE, provided_sex='male'),
            sample_row('CPG002', inferred_sex=FEMALE, provided_sex='female'),
            sample_row('CPG003', inferred_sex=MALE, provided_sex='female', gt_depth_mean=0.0),
        ],
        pair_rows=[
            pair_row('CPG001', 'CPG002', relatedness=0.01),
            pair_row('CPG001', 'CPG003', relatedness=0.02),
            pair_row('CPG002', 'CPG003', relatedness=0.03),
        ],
        ped=BROKEN_TRIO_PED,
    )

    flags, samples_df, _ = produce_flags(**inputs)

    # CPG003's sex mismatch is real but unusable, and no pair involving it may be flagged.
    assert 'CPG003' not in samples_df.sample_id.to_numpy()
    assert not [f for fs in flags.values() for f in fs if f.category == 'sex_inference_mismatch']
    assert all(
        'CPG003' not in (f.sg_id_1, f.sg_id_2)
        for fs in flags.values()
        for f in fs
        if f.category == 'relatedness_mismatch'
    )


def test_same_family_pairs_with_no_recorded_link_are_not_expected_unrelated(tmp_path):
    """
    Two family members with no blood path between them must not be flagged as 'expected unrelated'.

    CPG001 and CPG002 are both in FAM1 with no parents recorded, so peddy calls them 'unrelated'.
    The genotypes infer full siblings. Before the reframing this was a conflict claiming the
    pedigree said they were unrelated, which it never did.
    """
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            # somalier --infer reconstructs the sibship: both share the same two parents.
            sample_row('CPG001', inferred_sex=MALE, provided_sex='male', paternal='CPG004', maternal='CPG005'),
            sample_row('CPG002', inferred_sex=FEMALE, provided_sex='female', paternal='CPG004', maternal='CPG005'),
            sample_row('CPG004', inferred_sex=MALE, provided_sex='male'),
            sample_row('CPG005', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[
            pair_row('CPG001', 'CPG002', relatedness=0.49),
            pair_row('CPG001', 'CPG004', relatedness=0.48),
            pair_row('CPG001', 'CPG005', relatedness=0.51),
            pair_row('CPG002', 'CPG004', relatedness=0.5),
            pair_row('CPG002', 'CPG005', relatedness=0.49),
            pair_row('CPG004', 'CPG005', relatedness=0.01),
        ],
        # The expected pedigree knows only that all four are in FAM1, with no links at all.
        ped=(
            'FAM1\tCPG001\t0\t0\t1\t1\nFAM1\tCPG002\t0\t0\t2\t1\nFAM1\tCPG004\t0\t0\t1\t1\nFAM1\tCPG005\t0\t0\t2\t1\n'
        ),
    )

    flags, _, _ = produce_flags(**inputs)
    pedigree = [f for fs in flags.values() for f in fs if f.category == 'relatedness_mismatch']

    assert pedigree, 'expected the missing links to still be flagged, just not as "unrelated"'
    assert not [f for f in pedigree if f.expected_relationship == 'unrelated']
    assert all(f.expected_relationship == 'related at unknown level' for f in pedigree)


def test_cross_family_pairs_keep_expected_unrelated(tmp_path):
    """
    Two people in different families really are expected to be unrelated.

    So the reframing must not touch them: a related inference across a family boundary is the
    cross-family swap case, and it has to stay a conflict. Here somalier infers CPG001 and CPG002
    as full siblings while the expected pedigree has them in FAM1 and FAM2 respectively.
    """
    inputs = write_inputs(
        tmp_path,
        sample_rows=[
            sample_row('CPG001', inferred_sex=MALE, provided_sex='male', paternal='CPG004', maternal='CPG005'),
            sample_row('CPG002', inferred_sex=FEMALE, provided_sex='female', paternal='CPG004', maternal='CPG005'),
            sample_row('CPG004', inferred_sex=MALE, provided_sex='male'),
            sample_row('CPG005', inferred_sex=FEMALE, provided_sex='female'),
        ],
        pair_rows=[pair_row('CPG001', 'CPG002', relatedness=0.47)],
        ped=(
            'FAM1\tCPG001\t0\t0\t1\t1\nFAM2\tCPG002\t0\t0\t2\t1\nFAM1\tCPG004\t0\t0\t1\t1\nFAM1\tCPG005\t0\t0\t2\t1\n'
        ),
    )

    flags, _, _ = produce_flags(**inputs)
    pedigree = [f for fs in flags.values() for f in fs if f.category == 'relatedness_mismatch']

    assert [(f.expected_relationship, f.inferred_relationship) for f in pedigree] == [('unrelated', 'siblings')]
