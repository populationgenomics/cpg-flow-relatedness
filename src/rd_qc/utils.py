"""
suggested location for any utility methods or constants used across multiple stages
"""

from dataclasses import dataclass, field
from datetime import UTC, datetime
from functools import cache

from google.cloud import storage as gcs
from loguru import logger

from cpg_flow.metamist import get_metamist
from cpg_utils.config import config_retrieve, try_get_ar_guid
from metamist.graphql import gql, query


def get_gcs_object_size(fullpath: str, client: gcs.Client) -> int:
    """
    Get exact object size in GCS in GB, plus buffer for intermediate files.
    Returns 10 + buffer if the object is under 1GB.
    """
    buffer = config_retrieve(
        ['somalier_extract', 'storage_buffer'],
        20,
    )
    bucket_name, filepath = fullpath.removeprefix('gs://').split('/', 1)
    blob = client.bucket(bucket_name).blob(filepath)
    blob.reload()
    return max((blob.size // (1024**3), 10)) + buffer


SG_QUERY = gql("""
    query ProjectSomalier($project: String!) {
        project(name: $project) {
            sequencingGroups {
                id
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


@dataclass
class SgSomalierInfo:
    """Somalier fingerprint state for a single sequencing group."""

    sg_id: str
    participant_external_id: str
    somalier_path: str | None


@dataclass(kw_only=True)
class SomalierFlag:
    """Generic flag class for somalier QC checks."""

    category: str | None = None
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
    inferred_relationship: str
    relatedness: float
    ibs0: int
    ibs2: int


class SomalierIndex:
    """Dual-indexed view of somalier data: O(1) lookup by participant or sg_id."""

    def __init__(self, entries: list[SgSomalierInfo]):
        self.by_participant: dict[str, list[SgSomalierInfo]] = {}
        self.by_sg: dict[str, SgSomalierInfo] = {}
        for info in entries:
            self.by_participant.setdefault(info.participant_external_id, []).append(info)
            self.by_sg[info.sg_id] = info


@cache
def _query_project_sgs(project: str) -> list[dict]:
    """Cached metamist query — returns raw response data."""
    resolved = get_metamist().get_metamist_proj(project)
    response = query(SG_QUERY, variables={'project': resolved})
    return response['project']['sequencingGroups']


def get_project_sgs_and_fingerprints(project: str) -> SomalierIndex:
    """
    Query metamist for all SGs in the project with their somalier fingerprint status.
    Returns a SomalierIndex with O(1) lookup by participant or sg_id.

    Builds fresh SgSomalierInfo instances each call (safe to mutate)
    while the underlying metamist query is cached.
    """
    raw_sgs = _query_project_sgs(project)

    entries = []
    for sg in raw_sgs:
        sg_id = sg['id']
        participant_external_id = sg['sample']['participant']['externalId']
        analyses = sg.get('analyses', [])
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
    3. For unsequenced parents: keep external participant ID
    4. Substitute paternal_id/maternal_id with SG IDs where possible

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
