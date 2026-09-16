"""
Builders for synthetic somalier `relate` outputs, shared by every test that drives produce_flags.

The column lists mirror a real somalier 0.3.1 samples.tsv / pairs.tsv header exactly, because
`_check_sex` indexes into them by name and peddy re-reads samples.tsv as a PED file.

`sample_row` and `pair_row` are zero-argument-shaped builders in the sense the test guide means:
every field a test does not care about has a logically valid default, and a test overrides only
the fields its assertion depends on.
"""

from pathlib import Path

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

# somalier's own sex encoding in samples.tsv, and the sentinel it writes when inference failed.
MALE = 1
FEMALE = 2
SEX_INFERENCE_FAILED = 0
# What a PED carries when the provided sex was never recorded. `_check_sex` maps it to 'unknown'.
PROVIDED_SEX_UNKNOWN = '-9'

# somalier compares ~16.5k sites, and ibs0 is judged as a fraction of that, so tests must use a
# realistic count for the parent-child / siblings split to behave as it does in production.
SITES = 16500
IBS0_PARENT_CHILD = 0
IBS0_SIBLINGS = 340
IBS0_UNRELATED = 1414


def sample_row(
    sample_id: str,
    *,
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
    """The three `produce_flags` input paths, written under the test's own tmp_path."""
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


# TST001 (father) and TST002 (mother) are the recorded parents of TST003. peddy calls the two
# parents 'mom-dad', which expects unrelated, so a correct measurement flags nothing here.
# Columns are: family, sample, paternal, maternal, sex, phenotype.
TRIO_PED = 'FAM1\tTST001\t0\t0\t1\t1\nFAM1\tTST002\t0\t0\t2\t1\nFAM1\tTST003\tTST001\tTST002\t1\t1\n'


def write_pairs(tmp_path: Path, pair_rows: list[str], name: str = 'x.pairs.tsv') -> str:
    """Just a pairs.tsv, for callers that read it without the samples file or the PED."""
    pairs = tmp_path / name
    pairs.write_text('\t'.join(PAIR_COLUMNS) + '\n' + '\n'.join(pair_rows) + '\n')
    return str(pairs)
