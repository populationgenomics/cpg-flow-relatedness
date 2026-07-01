"""
suggested location for any utility methods or constants used across multiple stages
"""

from dataclasses import dataclass
import re
from functools import cache

from cpg_utils.config import config_retrieve
from loguru import logger
from metamist.graphql import gql, query


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
                    output
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
                analyses(type: {in_: ["cram", "gvcf", "vcf"]}) {
                    output
                    type
                    meta
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
    participant_id: str
    somalier_path: str | None

class SomalierIndex:
    """Dual-indexed view of somalier data: O(1) lookup by participant or sg_id."""

    def __init__(self, entries: list[SgSomalierInfo]):
        self.by_participant: dict[str, list[SgSomalierInfo]] = {}
        self.by_sg: dict[str, SgSomalierInfo] = {}
        for info in entries:
            self.by_participant.setdefault(info.participant_id, []).append(info)
            self.by_sg[info.sg_id] = info

def _resolve_project(project: str) -> str:
    if config_retrieve(['workflow', 'access_level']) == 'test' and not project.endswith('-test'):
        return project + '-test'
    return project

@cache
def get_project_sgs_and_fingerprints(project: str) -> dict[str, dict[str, str | None]]:
    """
    Query metamist for all SGs in the project, grouped by participant.

    Returns:
        {
            participant_external_id: {
                sg_id: somalier_filepath | None,
                ...
            },
            ...
        }
    """
    response = query(SG_QUERY, variables={'project': project})

    result: dict[str, dict[str, str | None]] = {}
    for sg in response['project']['sequencingGroups']:
        sg_id = sg['id']
        participant_id = sg['sample']['participant']['externalId']
        analyses = sg.get('analyses', [])
        somalier_path = analyses[0]['output'] if analyses else None

        if participant_id not in result:
            result[participant_id] = {}
        result[participant_id][sg_id] = somalier_path

    return result




def find_sgids_without_somalier(
    all_sgid_somaliers: dict[str, dict[str, str | None]],
) -> set[str]:
    """
    Find all SG IDs across all participants that don't have a somalier fingerprint.
    """
    missing = set()
    for sgid_map in all_sgid_somaliers.values():
        for sg_id, somalier_path in sgid_map.items():
            if somalier_path is None:
                missing.add(sg_id)
    return missing


def _select_best_file_for_sg(analyses: list[dict]) -> str | None:
    """
    Select the best source file for somalier extraction from a list of analyses.
    Default priority: CRAM > gVCF > VCF (configurable).
    """
    priority = config_retrieve(
        ['workflow', 'somalier_extract', 'priority'],
        ['cram', 'gvcf', 'vcf'],
    )

    buckets: dict[str, list[str]] = {'cram': [], 'gvcf': [], 'vcf': []}

    for analysis in analyses:
        output = analysis.get('output', '')
        if not output:
            continue
        if analysis.get('meta', {}).get('joint_called', False):
            continue

        analysis_type = analysis.get('type', '')
        if analysis_type == 'cram':
            buckets['cram'].append(output)
        elif analysis_type == 'gvcf':
            buckets['gvcf'].append(output)
        elif analysis_type == 'vcf':
            buckets['vcf'].append(output)

    for file_type in priority:
        if buckets.get(file_type):
            return buckets[file_type][0]
    return None

@cache
def select_somalier_extract_targets(project: str, sgids: tuple[str, ...]) -> dict[str, str]:
    """
    For each SG ID, query metamist for available analyses and select the best
    source file for somalier extraction.

    Returns {sg_id: source_file_path} for SGs where a suitable file was found.
    """
    resolved = _resolve_project(project)
    response = query(ANALYSIS_QUERY, variables={'project': resolved, 'sgIds': list(sgids)})

    targets: dict[str, str] = {}
    for sg in response['project']['sequencingGroups']:
        sg_id = sg['id']
        best_file = _select_best_file_for_sg(sg.get('analyses', []))
        if best_file:
            targets[sg_id] = best_file
        else:
            logger.warning(f'{sg_id}: no suitable file found for somalier extraction')

    return targets


@cache
def get_project_pedigree(project: str) -> list[dict]:
    """
    Query metamist for the full project pedigree.
    Returns the raw pedigree list with family_id, individual_id, paternal_id,
    maternal_id, sex, affected for every individual (including unsequenced).
    """
    resolved = _resolve_project(project)
    response = query(PEDIGREE_QUERY, variables={'project': resolved})
    return response['project']['pedigree']


def build_ped_content(
    project: str,
    sg_participant_map: dict[str, str],
) -> str:
    """
    Build a complete PED file string with SG IDs as individual identifiers.

    1. Query full pedigree from metamist (includes unsequenced parents)
    2. For participants with SGs: substitute SG ID for individual_id (one row per SG)
    3. For unsequenced parents: keep external participant ID
    4. Substitute paternal_id/maternal_id with SG IDs where possible

    Args:
        project: metamist project name
        sg_participant_map: {sg_id: external_participant_id}

    Returns:
        PED file content as a string (6-column format, tab-delimited)
    """
    pedigree = get_project_pedigree(project)

    # Build reverse mapping: participant_external_id -> [sg_id, ...]
    participant_to_sgs: dict[str, list[str]] = {}
    for sg_id, participant_id in sg_participant_map.items():
        participant_to_sgs.setdefault(participant_id, []).append(sg_id)

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

        sg_ids = participant_to_sgs.get(individual_id)
        if sg_ids:
            for sg_id in sorted(sg_ids):
                lines.append(f'{family_id}\t{sg_id}\t{paternal_id}\t{maternal_id}\t{sex}\t{affected}')
        else:
            lines.append(f'{family_id}\t{individual_id}\t{paternal_id}\t{maternal_id}\t{sex}\t{affected}')

    return '\n'.join(lines) + '\n'


def sg_ids_tag(sg_ids: list[str]) -> str:
    """Sorted, underscore-joined SG IDs for use in output file names."""
    return '_'.join(sorted(sg_ids))