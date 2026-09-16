"""
Mock Somalier flag data in the shape Metamist returns it, for offline report tests and renders.

Deliberately covers every branch the report's grouping logic has:

- a per-SG sex flag, active, and another resolved
- a sex flag whose provided sex is 'unknown'
- a pairwise self-relatedness flag
- a pedigree pair inside one family
- a pedigree pair straddling two families, which should render under both
- FAM02, which carries both active and resolved flags, so it appears in both sections
- TST009/TST011, which have no family, exercising the participant fallback
- TST001's sex flag, deliberately written without `sequencing_group_key`, so the report's
  derivation fallback for pre-existing flags gets exercised
- TST006 and TST007, which have no flags at all, so `total_sgs` exceeds the flagged count

`verdict` is spelled out per entry rather than computed with the real classifier: a fixture built
by the code under test can only show that the code round-trips its own output.

Nothing here touches Metamist. `MOCK_SEQUENCING_GROUPS` mirrors DATASET_SGS_QUERY's response and
`MOCK_SG_INFOS` mirrors what `get_sg_infos` would return, so a caller can skip both queries.
`test_get_sg_infos_*` drives the real parser against a recorded SGS_INFO_QUERY response, which is
what keeps this mirror honest.
"""

from rd_qc.scripts.somalier_flags_report import ReadFile, SGInfo
from rd_qc.utils import UNSPECIFIED_RELATED, VERDICT_CONFLICT, VERDICT_REFINEMENT

AR_GUID = 'mock-ar-guid-0001'
FIRST_SEEN = '2026-06-01T09:15:00+00:00'
RECENT = '2026-09-01T14:30:00+00:00'
RESOLVED_ON = '2026-08-22T11:05:00+00:00'
# The date part of FIRST_SEEN, which is what the read lines display.
READ_DATE = '2026-06-01'

SELF_RELATEDNESS_THRESHOLD = 0.9

# ---------------------------------------------------------------------------
# Identifiers the tests assert on by name
# ---------------------------------------------------------------------------
# FAM02 carries both active and resolved flags, and one end of the cross-family pair.
FAM_MIXED = 'FAM02'
# The other end of the cross-family pedigree pair.
FAM_CROSS = 'FAM07'
# Holds the self-relatedness pair and the unknown-provided-sex flag.
FAM_SELF = 'FAM03'
# Holds the one refinement, which belongs in the de-emphasised section.
FAM_REFINEMENT = 'FAM08'

PARTICIPANT_SELF = 'PID_B'
# Has no family at all, so its flags fall back to a participant-keyed group.
PARTICIPANT_NO_FAMILY = 'PID_G'

# The glance line leads with participants: TST004 is PID_C in FAM02, TST010 is PID_H in FAM07.
CROSS_FAMILY_SUBJECT = 'PID_C ↔ PID_H'
# Both members of this one are in FAM02.
SAME_FAMILY_SUBJECT = 'PID_C ↔ PID_D'

# Sizes and names the rendered read lines show.
CRAM_SIZE = '18.00 GiB'
FASTQ_SIZE = '1.50 GiB'
FIRST_SG_R1 = 'EXT_A_R1.fastq.gz'
# Recorded with neither a size nor a date, as some uploads are.
FIRST_SG_R2 = 'EXT_A_R2.fastq.gz'


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
    *,
    verdict: str,
    **overrides: object,
) -> dict:
    """
    A relatedness_mismatch flag as stored in SG meta.

    `inferred` is a measured degree from utils.infer_degree, not a peddy label. `verdict` is stated
    rather than derived, so this fixture is an independent oracle for the conflict/refinement
    split the report builds on — including the '' that pre-verdict records carry.
    """
    return {
        **_base('relatedness_mismatch', f'{sg_id_1}_{sg_id_2}'),
        'sg_id_1': sg_id_1,
        'sg_id_2': sg_id_2,
        'family_external_id': family,
        'expected_relationship': expected,
        'inferred_relationship': inferred,
        'verdict': verdict,
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
    'TST001': [legacy(sex_flag('TST001', provided='female', inferred='male', date=FIRST_SEEN))],
    # FAM03: an active self-relatedness pair for PID_B.
    'TST002': [self_relatedness_flag('TST002', 'TST003', participant=PARTICIPANT_SELF)],
    # FAM02: two active pedigree mismatches, one of them cross-family, plus one resolved.
    'TST004': [
        pedigree_flag(
            'TST004', 'TST005', FAM_MIXED, expected='parent-child', inferred='unrelated', verdict=VERDICT_CONFLICT
        ),
        pedigree_flag(
            'TST004',
            'TST010',
            FAM_MIXED,
            expected='unrelated',
            inferred='siblings',
            verdict=VERDICT_CONFLICT,
            relatedness=0.4911,
        ),
        resolved(
            pedigree_flag(
                'TST004', 'TST005', FAM_MIXED, expected='siblings', inferred='unrelated', verdict=VERDICT_CONFLICT
            )
        ),
    ],
    # FAM02 again: a resolved sex mismatch, so the family shows up in both sections.
    'TST005': [resolved(sex_flag('TST005', provided='male', inferred='female', date=FIRST_SEEN))],
    # FAM03: sex inference skipped upstream, so the provided sex is unknown.
    'TST008': [sex_flag('TST008', provided='unknown', inferred='male', x_sites=6, p_middling_ab=0.0712)],
    # No family at all: falls back to a participant-keyed group.
    'TST009': [
        self_relatedness_flag(
            'TST009', 'TST011', participant=PARTICIPANT_NO_FAMILY, relatedness=0.4218, ibs0=3902, ibs2=11244
        ),
    ],
    # FAM08: a refinement, not a conflict. The pedigree records no path between these two, and the
    # genotypes measure them as siblings, so it is a missing link rather than a contradiction.
    # Belongs in the de-emphasised section.
    'TST012': [
        pedigree_flag(
            'TST012',
            'TST013',
            FAM_REFINEMENT,
            expected=UNSPECIFIED_RELATED,
            inferred='siblings',
            verdict=VERDICT_REFINEMENT,
            relatedness=0.4873,
            ibs0=341,
        ),
    ],
}

_SG_IDENTITIES = [
    # sg_id, sample external id, participant, family, sample type
    ('TST001', 'EXT_A', 'PID_A', 'FAM01', 'blood'),
    ('TST002', 'EXT_B1', PARTICIPANT_SELF, FAM_SELF, 'blood'),
    ('TST003', 'EXT_B2', PARTICIPANT_SELF, FAM_SELF, 'saliva'),
    ('TST004', 'EXT_C', 'PID_C', FAM_MIXED, 'blood'),
    ('TST005', 'EXT_D', 'PID_D', FAM_MIXED, 'blood'),
    ('TST006', 'EXT_I', 'PID_I', 'FAM05', 'blood'),
    ('TST007', 'EXT_J', 'PID_J', 'FAM05', 'blood'),
    ('TST008', 'EXT_E', 'PID_E', FAM_SELF, 'blood'),
    ('TST009', 'EXT_G1', PARTICIPANT_NO_FAMILY, '', 'blood'),
    ('TST010', 'EXT_H', 'PID_H', FAM_CROSS, 'blood'),
    ('TST011', 'EXT_G2', PARTICIPANT_NO_FAMILY, '', 'saliva'),
    ('TST012', 'EXT_K', 'PID_K', FAM_REFINEMENT, 'blood'),
    ('TST013', 'EXT_L', 'PID_L', FAM_REFINEMENT, 'blood'),
]

# How many SGs the mock dataset holds, which is what `total_sgs` counts.
MOCK_TOTAL_SGS = len(_SG_IDENTITIES)

MOCK_SG_INFOS: dict[str, SGInfo] = {
    sg_id: SGInfo(
        sg_id=sg_id,
        sg_type='genome',
        sg_technology='short-read',
        sg_platform='illumina',
        crams=[ReadFile(f'{sample}.cram', size=CRAM_SIZE, date=READ_DATE)],
        fastq_pairs=[
            (
                ReadFile(f'{sample}_R1.fastq.gz', size=FASTQ_SIZE, date=READ_DATE),
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


# ---------------------------------------------------------------------------
# A response in the shape SGS_INFO_QUERY really returns
# ---------------------------------------------------------------------------
# Written in Metamist's own shape, including the quirks the parser exists to absorb: external IDs
# are a dict keyed by '' for the primary ID, `families` is a list even for one family and empty
# when there is none, a fastq pair arrives as one assay, and `datetime_added`/`size` can be null.
# `test_get_sg_infos_builds_the_same_infos_the_fixture_mirrors` drives the real parser over this
# and compares the result to MOCK_SG_INFOS, which is what stops the mirror drifting.
_CRAM_BYTES = 19327352832  # exactly 18 GiB
_FASTQ_R1_BYTES = 1610612736  # exactly 1.5 GiB


def _recorded_sg(sg_id: str, sample: str, participant: str, family: str, sample_type: str) -> dict:
    return {
        'id': sg_id,
        'meta': {},
        'type': 'genome',
        'technology': 'short-read',
        'platform': 'illumina',
        'assays': [
            {
                'id': f'assay-{sg_id}-cram',
                'meta': {
                    'reads_type': 'cram',
                    'reads': [
                        {
                            'basename': f'{sample}.cram',
                            'location': f'gs://test-bucket/{sample}.cram',
                            'size': _CRAM_BYTES,
                            'datetime_added': FIRST_SEEN,
                        },
                    ],
                },
            },
            {
                'id': f'assay-{sg_id}-fastq',
                'meta': {
                    'reads_type': 'fastq',
                    'reads': [
                        {
                            'basename': f'{sample}_R1.fastq.gz',
                            'location': f'gs://test-bucket/{sample}_R1.fastq.gz',
                            'size': _FASTQ_R1_BYTES,
                            'datetime_added': FIRST_SEEN,
                        },
                        {
                            'basename': f'{sample}_R2.fastq.gz',
                            'location': f'gs://test-bucket/{sample}_R2.fastq.gz',
                            'size': None,
                            'datetime_added': None,
                        },
                    ],
                },
            },
        ],
        'sample': {
            'id': f'XPGTEST{sg_id[-3:]}',
            'externalIds': {'': sample},
            'type': sample_type,
            'participant': {
                'id': 4100 + int(sg_id[-3:]),
                'externalIds': {'': participant},
                # Metamist omits the family entirely rather than nesting a blank one.
                'families': [{'id': 30 + int(sg_id[-3:]), 'externalIds': {'': family}}] if family else [],
            },
        },
    }


MOCK_SG_INFO_RESPONSE: dict = {
    'sequencingGroups': [_recorded_sg(*identity) for identity in _SG_IDENTITIES],
}
