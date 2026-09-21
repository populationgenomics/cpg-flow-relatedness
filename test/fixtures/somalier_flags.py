"""
Mock Somalier flag data in the shape Metamist returns it, for offline report tests and renders.

Deliberately covers every branch the report's grouping logic has:

- a per-SG sex flag, active, and another resolved
- a sex flag whose provided sex is 'unknown'
- a pairwise self-relatedness flag
- a pedigree pair inside one family
- a pedigree pair straddling two families, which should render under both
- FAM02, which carries both active and resolved flags, so it appears in both sections
- CPG009/CPG011, which have no family, exercising the participant fallback
- CPG001's sex flag, deliberately written without `sequencing_group_key`, so the report's
  derivation fallback for pre-existing flags gets exercised
- CPG006 and CPG007, which have no flags at all, so `total_sgs` exceeds the flagged count
- CPG014/CPG015, two sequencing groups of one participant, which the pedigree models as siblings

Nothing here touches Metamist. `MOCK_SEQUENCING_GROUPS` mirrors DATASET_SGS_QUERY's response and
`MOCK_SG_INFOS` mirrors what `get_sg_infos` would return, so a caller can skip both queries.
"""

from rd_qc.scripts.somalier_flags_report import ReadFile, SGInfo
from rd_qc.utils import UNSPECIFIED_RELATED, relatedness_verdict

AR_GUID = 'mock-ar-guid-0001'
FIRST_SEEN = '2026-06-01T09:15:00+00:00'
RECENT = '2026-09-01T14:30:00+00:00'
RESOLVED_ON = '2026-08-22T11:05:00+00:00'

SELF_RELATEDNESS_THRESHOLD = 0.9


def _base(category: str, sg_key: str, date: str = RECENT) -> dict:
    return {
        'category': category,
        'sequencing_group_key': sg_key,
        'date': date,
        'ar_guid': AR_GUID,
        'resolved': False,
        'resolution_date': None,
    }


def sex_flag(sg_key: str, provided: str, inferred: str, **overrides: object) -> dict:
    """A sex_inference_mismatch flag as stored in SG meta."""
    return {
        **_base('sex_inference_mismatch', sg_key),
        'provided': provided,
        'inferred': inferred,
        'mean_depth': 31.415926,
        'x_het_ratio': 0.016412,
        'x_depth_ratio': 1.0231,
        'y_depth_ratio': 0.98442,
        'x_sites': 4821,
        'p_middling_ab': 0.012377,
        **overrides,
    }


def self_relatedness_flag(sg_id_1: str, sg_id_2: str, participant: str, **overrides: object) -> dict:
    """A self_relatedness_mismatch flag as stored in SG meta."""
    return {
        **_base('self_relatedness_mismatch', f'{sg_id_1}_{sg_id_2}'),
        'sg_id_1': sg_id_1,
        'sg_id_2': sg_id_2,
        'participant_external_id': participant,
        'threshold': SELF_RELATEDNESS_THRESHOLD,
        'relatedness': 0.616722,
        'ibs0': 1204,
        'ibs2': 18337,
        **overrides,
    }


def pedigree_flag(
    sg_id_1: str,
    sg_id_2: str,
    family: str,
    expected: str,
    inferred: str,
    **overrides: object,
) -> dict:
    """
    A relatedness_mismatch flag as stored in SG meta.

    `inferred` is a measured degree from utils.infer_degree, not a peddy label. The verdict is
    computed with the real classifier rather than hardcoded, so a fixture can never disagree with
    the logic under test.
    """
    return {
        **_base('relatedness_mismatch', f'{sg_id_1}_{sg_id_2}'),
        'sg_id_1': sg_id_1,
        'sg_id_2': sg_id_2,
        'family_external_id': family,
        'expected_relationship': expected,
        'inferred_relationship': inferred,
        'verdict': relatedness_verdict(expected, inferred),
        'relatedness': 0.0214,
        'ibs0': 4821,
        'ibs2': 9033,
        **overrides,
    }


def resolved(flag: dict, resolution_date: str = RESOLVED_ON) -> dict:
    """Mark a flag resolved, the way reconciliation does."""
    return {**flag, 'resolved': True, 'resolution_date': resolution_date}


def legacy(flag: dict) -> dict:
    """Strip the SG key, as though the flag predates the sequencing_group_key field."""
    return {**flag, 'sequencing_group_key': ''}


# One entry per SG, keyed by SG id, holding only the flags that SG's meta carries. Pairwise flags
# are recorded against the first SG of the pair only, matching the producer scripts.
MOCK_FLAGS_BY_SG: dict[str, list[dict]] = {
    # FAM01: a single active sex mismatch, with no SG key so the fallback path runs.
    'CPG001': [legacy(sex_flag('CPG001', provided='female', inferred='male', date=FIRST_SEEN))],
    # FAM03: an active self-relatedness pair for PID_B.
    'CPG002': [self_relatedness_flag('CPG002', 'CPG003', participant='PID_B')],
    # FAM02: two active pedigree mismatches, one of them cross-family, plus one resolved.
    'CPG004': [
        pedigree_flag('CPG004', 'CPG005', 'FAM02', expected='parent-child', inferred='unrelated'),
        pedigree_flag('CPG004', 'CPG010', 'FAM02', expected='unrelated', inferred='siblings', relatedness=0.4911),
        resolved(pedigree_flag('CPG004', 'CPG005', 'FAM02', expected='siblings', inferred='unrelated')),
    ],
    # FAM02 again: a resolved sex mismatch, so the family shows up in both sections.
    'CPG005': [resolved(sex_flag('CPG005', provided='male', inferred='female', date=FIRST_SEEN))],
    # FAM03: sex inference skipped upstream, so the provided sex is unknown.
    'CPG008': [sex_flag('CPG008', provided='unknown', inferred='male', x_sites=6, p_middling_ab=0.0712)],
    # No family at all: falls back to a participant-keyed group.
    'CPG009': [
        self_relatedness_flag('CPG009', 'CPG011', participant='PID_G', relatedness=0.4218, ibs0=3902, ibs2=11244),
    ],
    # FAM08: a refinement, not a conflict. The pedigree records no path between these two, and the
    # genotypes measure them as siblings, so it is a missing link rather than a contradiction.
    # Belongs in the de-emphasised section.
    'CPG012': [
        pedigree_flag(
            'CPG012',
            'CPG013',
            'FAM08',
            expected=UNSPECIFIED_RELATED,
            inferred='siblings',
            relatedness=0.4873,
            ibs0=341,
        ),
    ],
    # FAM09: two sequencing groups for one person (PID_M). build_ped_content writes a PED row per
    # SG, both with the same parents, so peddy calls them full siblings while the genotypes
    # measure them as identical. Mirrors KDG0134PR in ghfm-kidgen.
    'CPG014': [
        pedigree_flag(
            'CPG014',
            'CPG015',
            'FAM09',
            expected='full siblings',
            inferred='identical',
            relatedness=0.9982,
            ibs0=3,
            ibs2=19871,
        ),
    ],
}

_SG_IDENTITIES = [
    # sg_id, sample external id, participant, family, sample type
    ('CPG001', 'EXT_A', 'PID_A', 'FAM01', 'blood'),
    ('CPG002', 'EXT_B1', 'PID_B', 'FAM03', 'blood'),
    ('CPG003', 'EXT_B2', 'PID_B', 'FAM03', 'saliva'),
    ('CPG004', 'EXT_C', 'PID_C', 'FAM02', 'blood'),
    ('CPG005', 'EXT_D', 'PID_D', 'FAM02', 'blood'),
    ('CPG006', 'EXT_I', 'PID_I', 'FAM05', 'blood'),
    ('CPG007', 'EXT_J', 'PID_J', 'FAM05', 'blood'),
    ('CPG008', 'EXT_E', 'PID_E', 'FAM03', 'blood'),
    ('CPG009', 'EXT_G1', 'PID_G', '', 'blood'),
    ('CPG010', 'EXT_H', 'PID_H', 'FAM07', 'blood'),
    ('CPG011', 'EXT_G2', 'PID_G', '', 'saliva'),
    ('CPG012', 'EXT_K', 'PID_K', 'FAM08', 'blood'),
    ('CPG013', 'EXT_L', 'PID_L', 'FAM08', 'blood'),
    # Two SGs, one participant: the same-individual case.
    ('CPG014', 'EXT_M1', 'PID_M', 'FAM09', 'blood'),
    ('CPG015', 'EXT_M2', 'PID_M', 'FAM09', 'saliva'),
]

MOCK_SG_INFOS: dict[str, SGInfo] = {
    sg_id: SGInfo(
        sg_id=sg_id,
        sg_type='genome',
        sg_technology='short-read',
        sg_platform='illumina',
        crams=[ReadFile(f'{sample}.cram', size='18.00 GiB', date=FIRST_SEEN[:10])],
        fastq_pairs=[
            (
                ReadFile(f'{sample}_R1.fastq.gz', size='1.50 GiB', date=FIRST_SEEN[:10]),
                # No size or date, as some uploads record neither.
                ReadFile(f'{sample}_R2.fastq.gz'),
            ),
        ],
        other_reads=[],
        sample_external_id=sample,
        sample_type=sample_type,
        participant_external_id=participant,
        family_external_id=family,
    )
    for sg_id, sample, participant, family, sample_type in _SG_IDENTITIES
}

MOCK_SEQUENCING_GROUPS: list[dict] = [
    {'id': sg_id, 'type': 'genome', 'meta': {'somalier_flags': MOCK_FLAGS_BY_SG.get(sg_id, [])}}
    for sg_id, *_ in _SG_IDENTITIES
]

# A second dataset with no flags at all, for the all-clear render.
MOCK_ALL_CLEAR_SEQUENCING_GROUPS: list[dict] = [
    {'id': sg_id, 'type': 'genome', 'meta': {}} for sg_id, *_ in _SG_IDENTITIES
]
