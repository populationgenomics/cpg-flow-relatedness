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
                        id
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

SOMALIER_RELATE_ANALYSES_QUERY = gql("""
    query SgSomalierRelateAnalyses($project: String!) {
        project(name: $project) {
            analyses(type: {in_: ["somalier_relate"]}) {
                outputs
                type
                meta
                sequencingGroups {
                    id
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
    participant_id: int
    participant_external_id: str
    somalier_path: str | Path | None


class SomalierIndex:
    """Dual-indexed view of somalier data: O(1) lookup by participant or sg_id."""

    def __init__(self, entries: list[SgSomalierInfo]):
        self.by_participant: dict[tuple[int, str], list[SgSomalierInfo]] = {}
        self.by_sg: dict[str, SgSomalierInfo] = {}
        for info in entries:
            self.by_participant.setdefault((info.participant_id, info.participant_external_id), []).append(info)
            self.by_sg[info.sg_id] = info


@dataclass(kw_only=True)
class SomalierFlag:
    """Generic flag class for somalier QC checks."""

    category: str | None = None
    # Sorted, underscore-joined SG IDs this flag involves: 'CPG001' for per-SG flags,
    # 'CPG002_CPG003' for pairwise ones. Set during reconciliation in record_somalier_flags.py.
    # Defaults to '' so older records without the field still deserialise; readers fall back to
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
    # 'conflict' or 'refinement', from relatedness_verdict. Stored so the report can classify a
    # flag without the pedigree. Empty on older records; readers treat those as conflicts.
    verdict: str = ''


# peddy's Ped.relation() vocabulary. Pinned here because the pedigree check compares these exact
# strings, so an upstream rename would silently reclassify flags.
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

# What peddy reports for a pair in the same family with no recorded path between them. Distinct
# from 'unknown', which is our own fallback for a sample missing from the PED entirely.
UNSPECIFIED_RELATED = 'related at unknown level'


def refine_expected_relationship(relation: str, family_1: str | None, family_2: str | None) -> str:
    """
    Reinterpret peddy's 'unrelated' for two individuals recorded in the same family.

    peddy returns 'unrelated' for any pair with no blood path, including two members of one family
    whose connecting links are simply not recorded. For those, the honest expectation is
    relatedness at some unspecified level rather than "expected unrelated".

    Only the name of the expectation changes; `relatedness_verdict` decides whether it is met.
    Co-parents are unaffected, as peddy reports a recorded mother and father as 'mom-dad'.
    """
    if relation == 'unrelated' and family_1 and family_1 == family_2:
        return UNSPECIFIED_RELATED
    return relation


# Relatedness degrees, derived from the kinship coefficient rather than from somalier's `--infer`
# pedigree reconstruction. --infer assumes high quality pairs with both parents present; pedigrees
# with only one recorded parent make it rewrite family IDs and invent parent links, which any
# downstream pedigree check then reads as truth.
DEGREE_IDENTICAL = 'identical'
DEGREE_PARENT_CHILD = 'parent-child'
DEGREE_SIBLINGS = 'siblings'
DEGREE_SECOND = 'second-degree'
DEGREE_THIRD = 'third-degree'
DEGREE_UNRELATED = 'unrelated'

# Kinship coefficient lower bounds. somalier reports relatedness on the 2*phi scale, where each
# successive degree halves: identical 1.0, first-degree 0.5, second 0.25, third 0.125. The bounds
# are the geometric midpoints between those expectations.
IDENTICAL_MIN_RELATEDNESS = 0.90
FIRST_DEGREE_MIN_RELATEDNESS = 0.354
SECOND_DEGREE_MIN_RELATEDNESS = 0.177
THIRD_DEGREE_MIN_RELATEDNESS = 0.088

# ibs0 splits the two first-degree relationships: a parent and child share an allele at every site,
# so ibs0 is ~0, while full siblings inherit different alleles at some. Expressed as a fraction of
# the sites compared so the threshold survives a different sites VCF.
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


# Measured degrees consistent with each relationship a pedigree can state. peddy's 'siblings' means
# "shares at least one recorded parent", so it spans full and half siblings.
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
    unspecified expectation means the pedigree records no path between the pair. Any relatedness
    there is a plausible missing link, so it is a refinement -- except two identical genomes, which
    indicate a duplicate or swap rather than an omission. Measuring unrelated against an
    unspecified expectation says nothing, so it is not flagged.
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
    # Third-degree relatedness against an expectation of unrelated is not assertable: it overlaps
    # the cohort background distribution, so surface it as a refinement rather than an error.
    if measured == DEGREE_THIRD and acceptable == frozenset({DEGREE_UNRELATED}):
        return VERDICT_REFINEMENT
    return VERDICT_CONFLICT


def convert_to_web_url(dataset_name: str, html_path: Path | str) -> str:
    """Convert a gs:// web-bucket path to the http(s) web URL."""
    # storage config is keyed on the main dataset, not the -test suffix
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

    The metamist query is cached, but fresh SgSomalierInfo instances are built each call.
    With filter_sgs, only SGs matching the configured sequencing type and technology are included.
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
        participant_id = sg['sample']['participant']['id']
        participant_external_id = sg['sample']['participant']['externalId']
        analyses = sg['analyses']
        somalier_path = analyses[0]['outputs'].get('path') if analyses else None
        entries.append(SgSomalierInfo(sg_id, participant_id, participant_external_id, somalier_path))

    return SomalierIndex(entries)


def find_sgids_without_somalier(index: SomalierIndex) -> set[str]:
    """Find all SG IDs that don't have a somalier fingerprint."""
    return {info.sg_id for info in index.by_sg.values() if info.somalier_path is None}


def _select_best_file_for_sg(analyses: list[dict]) -> str | None:
    """Pick the highest-priority analysis output, CRAM before gVCF by default."""
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
    Select the best somalier extraction source file for each SG ID.

    Returns ({sg_id: source_file_path}, {sg_id: external IDs}) for SGs with a suitable file.
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
def get_somalier_relate_analyses(project: str) -> dict[str, list[dict]]:
    """
    Retrieve all somalier relate analyses for the given project.
    """
    resolved = get_metamist().get_metamist_proj(project)
    response = query(SOMALIER_RELATE_ANALYSES_QUERY, variables={'project': resolved})
    analyses_by_sg: dict[str, list[dict]] = {}
    for analysis in response['project']['analyses']:
        sgs = analysis['sequencingGroups']
        for sg in sgs:
            sg_id = sg['id']
            analyses_by_sg.setdefault(sg_id, []).append(analysis)
    return analyses_by_sg

@cache
def get_project_pedigree(project: str) -> list[dict]:
    """
    Query metamist for the full project pedigree, including unsequenced individuals.

    Each row has family_id, individual_id, paternal_id, maternal_id, sex and affected.
    """
    resolved = get_metamist().get_metamist_proj(project)
    response = query(PEDIGREE_QUERY, variables={'project': resolved})
    return response['project']['pedigree']


def build_ped_content(
    project: str,
    index: SomalierIndex,
) -> str:
    """
    Build tab-delimited 6-column PED content, using SG IDs as individual identifiers.

    Participants with SGs get one row per SG; unsequenced individuals keep their participant ID.
    Paternal/maternal IDs are substituted with SG IDs where possible.
    """
    pedigree = get_project_pedigree(project)

    participant_to_sgs: dict[str, list[str]] = {}
    for info in index.by_sg.values():
        participant_to_sgs.setdefault(info.participant_external_id, []).append(info.sg_id)

    # parent columns hold a single ID, so pick one SG per participant
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
