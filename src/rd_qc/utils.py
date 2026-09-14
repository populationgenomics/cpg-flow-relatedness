"""
suggested location for any utility methods or constants used across multiple stages
"""

from dataclasses import dataclass, field
from datetime import UTC, datetime
from functools import cache

from google.cloud import storage as gcs
from loguru import logger

from cpg_flow.metamist import get_metamist
from cpg_utils import Path
from cpg_utils.config import config_retrieve, try_get_ar_guid
from metamist.graphql import gql, query

GCS_CLIENT: gcs.Client | None = None


SG_QUERY = gql("""
    query ProjectSomalier($project: String!) {
        project(name: $project) {
            sequencingGroups {
                id
                type
                technology
                sample {
                    participant {
                        externalId
                    }
                }
                analyses(type: {eq: "somalier"}) {
                    outputs
                    meta
                }
            }
        }
    }
""")

ANALYSIS_QUERY = gql("""
    query SgAnalyses($project: String!, $sgIds: [String!]!) {
        project(name: $project) {
            sequencingGroups(id: {in_: $sgIds}) {
                id
                analyses(type: {in_: ["cram", "gvcf"]}) {
                    outputs
                    type
                    meta
                }
                sample {
                    externalId
                    participant {
                        externalId
                    }
                }
            }
        }
    }
""")

PEDIGREE_QUERY = gql("""
    query ProjectPedigree($project: String!) {
        project(name: $project) {
            pedigree
        }
    }
""")


def get_gcs_client():
    global GCS_CLIENT
    if GCS_CLIENT is None:
        GCS_CLIENT = gcs.Client()
    return GCS_CLIENT


def get_gcs_object_size(fullpath: str) -> int:
    """
    Get exact object size in GCS in GB, plus buffer for intermediate files.
    Returns 10 + buffer if the object is under 1GB.
    """
    buffer = config_retrieve(
        ['workflow', 'somalier_extract', 'storage_buffer'],
        20,
    )
    bucket_name, filepath = fullpath.removeprefix('gs://').split('/', 1)
    blob = get_gcs_client().bucket(bucket_name).blob(filepath)
    blob.reload()
    return max((blob.size // (1024**3), 10)) + buffer


@dataclass
class SgSomalierInfo:
    """Somalier fingerprint state for a single sequencing group."""

    sg_id: str
    participant_external_id: str
    somalier_path: str | Path | None


class SomalierIndex:
    """Dual-indexed view of somalier data: O(1) lookup by participant or sg_id."""

    def __init__(self, entries: list[SgSomalierInfo]):
        self.by_participant: dict[str, list[SgSomalierInfo]] = {}
        self.by_sg: dict[str, SgSomalierInfo] = {}
        for info in entries:
            self.by_participant.setdefault(info.participant_external_id, []).append(info)
            self.by_sg[info.sg_id] = info


@dataclass(kw_only=True)
class SomalierFlag:
    """Generic flag class for somalier QC checks."""

    category: str | None = None
    # Sorted, underscore-joined SG IDs this flag involves: 'CPG001' for per-SG flags,
    # 'CPG002_CPG003' for pairwise ones. Set during reconciliation in record_somalier_flags.py,
    # which is also where the per-category identity keys live. Defaults to '' so that flags
    # already recorded in Metamist without the field still deserialise; readers fall back to
    # deriving it from sg_id_1/sg_id_2.
    sequencing_group_key: str = ''
    date: str = field(default_factory=lambda: datetime.now(tz=UTC).isoformat(timespec='seconds'))
    ar_guid: str = field(default_factory=try_get_ar_guid)
    resolved: bool = False
    resolution_date: str | None = None


@dataclass(kw_only=True)
class SomalierSexInferenceFlag(SomalierFlag):
    """Somalier sex inference mismatch flag."""

    provided: str
    inferred: str
    mean_depth: float
    x_het_ratio: float  # the actual decision statistic
    x_depth_ratio: float  # ~1 XY, ~2 XX
    y_depth_ratio: float  # ~1 XY, ~0 YY
    x_sites: int  # <10 => no call attempted
    p_middling_ab: float  # >=0.06 => inference skipped


@dataclass(kw_only=True)
class SomalierSelfRelatednessFlag(SomalierFlag):
    """Somalier self-relatedness mismatch flag."""

    sg_id_1: str
    sg_id_2: str
    participant_external_id: str
    threshold: float
    relatedness: float
    ibs0: int
    ibs2: int


@dataclass(kw_only=True)
class SomalierRelatednessFlag(SomalierFlag):
    """Somalier relatedness mismatch flag."""

    sg_id_1: str
    sg_id_2: str
    family_external_id: str
    expected_relationship: str
    # The degree the measurement supports, from infer_degree. Not somalier's --infer output.
    inferred_relationship: str
    relatedness: float
    ibs0: int
    ibs2: int
    # 'conflict' or 'refinement', from relatedness_verdict. Stored rather than re-derived so the
    # report does not need the pedigree to classify a flag. Empty on flags recorded before this
    # field existed; readers fall back to treating those as conflicts.
    verdict: str = ''


# peddy's Ped.relation() vocabulary. Pinned here because the pedigree check classifies flags by
# comparing these exact strings, so an upstream rename would silently reclassify everything.
PEDDY_RELATIONSHIPS = frozenset(
    {
        'cousins',
        'full siblings',
        'grandchild',
        'mom-dad',
        'niece/nephew',
        'parent-child',
        'related at unknown level',
        'siblings',
        'unknown',
        'unrelated',
    }
)

# What peddy reports when the pedigree puts a pair in the same family but records no path between
# them. Note it is not the string 'unknown', which is what our own code falls back to when a sample
# is missing from the PED entirely.
UNSPECIFIED_RELATED = 'related at unknown level'


def refine_expected_relationship(relation: str, family_1: str | None, family_2: str | None) -> str:
    """
    Reinterpret peddy's 'unrelated' for two individuals recorded in the same family.

    peddy returns 'unrelated' whenever it finds no blood path between a pair, which includes two
    members of one family whose connecting links simply are not recorded. Reporting that as
    "expected unrelated" is wrong: all 66 of perth-neuro's 'expected unrelated' flags were
    same-family and not one was cross-family. If the pedigree puts two people in a family, the
    honest expectation is that they are related at some unspecified level.

    This only changes what the expectation is *called*. Whether that expectation is met is decided
    by `relatedness_verdict`, which treats an unspecified expectation as satisfied by anything
    short of two identical samples. So a pair of in-laws whom the genotypes also call unrelated
    stays unflagged.

    Co-parents are unaffected, because peddy reports a recorded mother and father as 'mom-dad'
    rather than 'unrelated'. A mother and father who turn out to be blood relatives therefore
    still surfaces as a conflict rather than being quietly demoted.
    """
    if relation == 'unrelated' and family_1 and family_1 == family_2:
        return UNSPECIFIED_RELATED
    return relation


# ---------------------------------------------------------------------------
# Relatedness degrees, inferred from what somalier measured
# ---------------------------------------------------------------------------
# We derive the relationship ourselves from the kinship coefficient rather than reading somalier's
# `--infer` reconstruction. somalier documents --infer as being for high quality sample pairs where
# both parents are present, and CPG pedigrees frequently record only one parent, so it runs well
# outside its envelope: on perth-neuro it renumbered the family id of 477 of 623 samples, invented
# 114 parent links, and created 42 synthetic placeholder parents. peddy then faithfully read that
# fabricated pedigree, which produced ~160 false flags while missing a genotypically identical pair.
DEGREE_IDENTICAL = 'identical'
DEGREE_PARENT_CHILD = 'parent-child'
DEGREE_SIBLINGS = 'siblings'
DEGREE_SECOND = 'second-degree'
DEGREE_THIRD = 'third-degree'
DEGREE_UNRELATED = 'unrelated'

# Kinship coefficient lower bounds. somalier reports relatedness on the 2*phi scale, where each
# successive degree halves: identical 1.0, first-degree 0.5, second 0.25, third 0.125. The bounds
# are the geometric midpoints between those expectations, which is both the principled split and
# what perth-neuro's measured clusters support: parent-child 0.438..0.548 (n=215), full siblings
# 0.430..0.554 (n=26), grandchild 0.203..0.306 (n=8), niece/nephew 0.204..0.296 (n=8).
#
# The second-degree bound matters most. Cross-family pairs, which the pedigree really does expect
# to be unrelated, form a smooth background distribution with no upper cluster: median -0.006,
# p99.9 0.072, max 0.158 over 99,515 pairs. So 0.177 sits in the genuine gap between that
# background and the real second-degree cluster at 0.203+.
IDENTICAL_MIN_RELATEDNESS = 0.90
FIRST_DEGREE_MIN_RELATEDNESS = 0.354
SECOND_DEGREE_MIN_RELATEDNESS = 0.177
THIRD_DEGREE_MIN_RELATEDNESS = 0.088

# ibs0 splits the two first-degree relationships: a parent and child share an allele at every site,
# so ibs0 is ~0, while full siblings inherit different alleles at some. Expressed as a fraction of
# the sites compared so the threshold survives a different sites VCF. On perth-neuro parent-child
# reached at most 0.0007 (ibs0 <= 12) and siblings never went below 0.012 (ibs0 >= 201), a clean
# separation with no overlap, so this sits between them with margin on both sides.
PARENT_CHILD_MAX_IBS0_RATIO = 0.005


def infer_degree(relatedness: float, ibs0: int, sites: int) -> str:
    """
    The relatedness degree the measurement supports, independent of any pedigree.

    `sites` is somalier's `n`, the number of sites the pair was compared at.
    """
    if relatedness >= IDENTICAL_MIN_RELATEDNESS:
        return DEGREE_IDENTICAL
    if relatedness >= FIRST_DEGREE_MIN_RELATEDNESS:
        ratio = (ibs0 / sites) if sites else 0.0
        return DEGREE_PARENT_CHILD if ratio <= PARENT_CHILD_MAX_IBS0_RATIO else DEGREE_SIBLINGS
    if relatedness >= SECOND_DEGREE_MIN_RELATEDNESS:
        return DEGREE_SECOND
    if relatedness >= THIRD_DEGREE_MIN_RELATEDNESS:
        return DEGREE_THIRD
    return DEGREE_UNRELATED


# Which measured degrees are consistent with each relationship a pedigree can state. peddy's
# 'siblings' means "shares at least one recorded parent", so it spans full and half siblings and
# accepts either a first- or second-degree measurement.
EXPECTED_DEGREES: dict[str, frozenset[str]] = {
    'parent-child': frozenset({DEGREE_PARENT_CHILD}),
    'full siblings': frozenset({DEGREE_SIBLINGS}),
    'siblings': frozenset({DEGREE_SIBLINGS, DEGREE_SECOND}),
    'grandchild': frozenset({DEGREE_SECOND}),
    'niece/nephew': frozenset({DEGREE_SECOND}),
    'cousins': frozenset({DEGREE_THIRD}),
    'mom-dad': frozenset({DEGREE_UNRELATED}),
    'unrelated': frozenset({DEGREE_UNRELATED}),
}

VERDICT_OK = 'ok'
VERDICT_CONFLICT = 'conflict'
VERDICT_REFINEMENT = 'refinement'


def expected_degrees(relationship: str) -> frozenset[str] | None:
    """
    The measured degrees consistent with the stated relationship, or None when the pedigree states
    nothing usable ('related at unknown level', or a sample missing from the pedigree entirely).
    """
    return EXPECTED_DEGREES.get(relationship)


def relatedness_verdict(relationship: str, measured: str) -> str:
    """
    Whether a measured degree agrees with the pedigree, contradicts it, or fills a gap in it.

    `relationship` is the expected relationship after `refine_expected_relationship`, so an
    unspecified expectation means the pedigree records no path between the pair. In that case any
    degree of relatedness is a plausible missing link and counts as a refinement, with one
    exception: two identical genomes are never a pedigree omission, they are one sample recorded
    twice or a swap, so that stays a conflict. Measuring unrelated against an unspecified
    expectation says nothing at all, which is the in-law case, so it is not flagged.
    """
    acceptable = expected_degrees(relationship)
    if acceptable is None:
        if measured == DEGREE_IDENTICAL:
            return VERDICT_CONFLICT
        if measured == DEGREE_UNRELATED:
            return VERDICT_OK
        return VERDICT_REFINEMENT
    if measured in acceptable:
        return VERDICT_OK
    # A third-degree measurement against an expectation of unrelated is not assertable. Distant
    # relatedness is indistinguishable from cohort background: on perth-neuro the cross-family
    # background reached 0.158 and its p99.99 was 0.126, which is exactly where a real first
    # cousin sits. The 40 pairs this catches were also concentrated on a handful of samples, one
    # of them appearing in 8 different pairs, which is the signature of a sample-level artefact
    # rather than kinship. So it is surfaced as a refinement rather than claimed as an error.
    if measured == DEGREE_THIRD and acceptable == frozenset({DEGREE_UNRELATED}):
        return VERDICT_REFINEMENT
    return VERDICT_CONFLICT


def convert_to_web_url(dataset_name: str, html_path: Path | str) -> str:
    """
    Convert a gs:// web-bucket path to the http(s) web URL.
    """
    # Important - strip -test from dataset suffix before constructing the web URL
    dataset_name = dataset_name.removesuffix('-test')
    return str(html_path).replace(
        config_retrieve(['storage', dataset_name, 'web']),
        config_retrieve(['storage', dataset_name, 'web_url']),
    )


@cache
def _query_project_sgs(project: str) -> list[dict]:
    """Cached metamist query — returns raw response data."""
    resolved = get_metamist().get_metamist_proj(project)
    response = query(SG_QUERY, variables={'project': resolved})
    return response['project']['sequencingGroups']


def get_project_sgs_and_fingerprints(project: str, filter_sgs: bool = False) -> SomalierIndex:
    """
    Query metamist for all SGs in the project with their somalier fingerprint status.
    Returns a SomalierIndex with O(1) lookup by participant or sg_id.

    Builds fresh SgSomalierInfo instances each call (safe to mutate)
    while the underlying metamist query is cached.

    If filter_sgs is True, only include SGs that meet the sequencing type & technology requirements
    as defined in the config.
    """
    raw_sgs = _query_project_sgs(project)

    entries = []
    for sg in raw_sgs:
        if filter_sgs:
            seq_type = sg['type']
            seq_tech = sg['technology']
            if seq_type != config_retrieve(['workflow', 'sequencing_type']):
                logger.debug(f'{sg["id"]}: skipping SG with sequencing type {seq_type}')
                continue
            if seq_tech != config_retrieve(['workflow', 'sequencing_technology']):
                logger.debug(f'{sg["id"]}: skipping SG with sequencing technology {seq_tech}')
                continue
        sg_id = sg['id']
        participant_external_id = sg['sample']['participant']['externalId']
        analyses = sg['analyses']
        somalier_path = analyses[0]['outputs'].get('path') if analyses else None
        entries.append(SgSomalierInfo(sg_id, participant_external_id, somalier_path))

    return SomalierIndex(entries)


def find_sgids_without_somalier(index: SomalierIndex) -> set[str]:
    """Find all SG IDs that don't have a somalier fingerprint."""
    return {info.sg_id for info in index.by_sg.values() if info.somalier_path is None}


def _select_best_file_for_sg(analyses: list[dict]) -> str | None:
    """
    Select the best source file for somalier extraction from a list of analyses.
    Default priority: CRAM > gVCF > VCF (configurable).
    """
    priority = config_retrieve(
        ['somalier_extract', 'priority'],
        ['cram', 'gvcf'],
    )

    buckets: dict[str, list[str]] = {'cram': [], 'gvcf': []}

    for analysis in analyses:
        output = (analysis.get('outputs') or {}).get('path', '')
        if not output:
            continue
        if analysis.get('meta', {}).get('joint_called', False):
            continue

        analysis_type = analysis.get('type', '')
        if analysis_type == 'cram':
            buckets['cram'].append(output)
        elif analysis_type == 'gvcf':
            buckets['gvcf'].append(output)

    for file_type in priority:
        if buckets.get(file_type):
            return buckets[file_type][0]
    return None


@cache
def select_somalier_extract_targets(project: str, sgids: tuple[str, ...]) -> tuple[dict[str, str], dict[str, dict]]:
    """
    For each SG ID, query metamist for available analyses and select the best
    source file for somalier extraction.

    Returns {sg_id: source_file_path} for SGs where a suitable file was found.
    """
    resolved = get_metamist().get_metamist_proj(project)
    response = query(ANALYSIS_QUERY, variables={'project': resolved, 'sgIds': list(sgids)})

    sg_id_map: dict[str, dict] = {}
    targets: dict[str, str] = {}
    for sg in response['project']['sequencingGroups']:
        sg_id = sg['id']
        best_file = _select_best_file_for_sg(sg.get('analyses', []))
        if best_file:
            targets[sg_id] = best_file
            sg_id_map[sg_id] = {
                'sample_external_id': sg['sample']['externalId'],
                'participant_external_id': sg['sample']['participant']['externalId'],
            }
        else:
            logger.warning(f'{sg_id}: no suitable file found for somalier extraction')

    return targets, sg_id_map


@cache
def get_project_pedigree(project: str) -> list[dict]:
    """
    Query metamist for the full project pedigree.
    Returns the raw pedigree list with family_id, individual_id, paternal_id,
    maternal_id, sex, affected for every individual (including unsequenced).
    """
    resolved = get_metamist().get_metamist_proj(project)
    response = query(PEDIGREE_QUERY, variables={'project': resolved})
    return response['project']['pedigree']


def build_ped_content(
    project: str,
    index: SomalierIndex,
) -> str:
    """
    Build a complete PED file string with SG IDs as individual identifiers.

    1. Query full pedigree from metamist (includes unsequenced parents)
    2. For participants with SGs: substitute SG ID for individual_id (one row per SG)
    3. Substitute paternal_id/maternal_id with SG IDs where possible, otherwise leave as participant ID

    Args:
        project: metamist project name
        index: SomalierIndex with participant and SG data

    Returns:
        PED file content as a string (6-column format, tab-delimited)
    """
    pedigree = get_project_pedigree(project)

    # Build reverse mapping: participant_external_id -> [sg_id, ...]
    participant_to_sgs: dict[str, list[str]] = {}
    for info in index.by_sg.values():
        participant_to_sgs.setdefault(info.participant_external_id, []).append(info.sg_id)

    # For ID substitution in paternal/maternal fields, pick first SG per participant
    participant_to_primary_sg: dict[str, str] = {}
    for participant_id, sg_ids in participant_to_sgs.items():
        participant_to_primary_sg[participant_id] = sorted(sg_ids)[0]

    def _resolve_id(individual_id: str | None) -> str | None:
        if individual_id is None:
            return None
        return participant_to_primary_sg.get(individual_id, individual_id)

    lines = []
    for row in pedigree:
        individual_id = row['individual_id']
        family_id = row['family_id']
        paternal_id = _resolve_id(row.get('paternal_id')) or '0'
        maternal_id = _resolve_id(row.get('maternal_id')) or '0'
        sex = row.get('sex', 0)
        affected = row.get('affected', -9)

        sg_id_list = participant_to_sgs.get(individual_id)
        if sg_id_list:
            for sg_id in sorted(sg_id_list):
                lines.append(f'{family_id}\t{sg_id}\t{paternal_id}\t{maternal_id}\t{sex}\t{affected}')
        else:
            lines.append(f'{family_id}\t{individual_id}\t{paternal_id}\t{maternal_id}\t{sex}\t{affected}')

    return '\n'.join(lines) + '\n'


def sg_ids_tag(sg_ids: list[str]) -> str:
    """Sorted, underscore-joined SG IDs for use in output file names."""
    return '_'.join(sorted(sg_ids))
