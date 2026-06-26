"""
suggested location for any utility methods or constants used across multiple stages
"""

from dataclasses import dataclass
import re
from functools import cache

from cpg_utils.config import config_retrieve
from loguru import logger
from metamist.graphql import gql, query

_CRAM_PATTERN = re.compile(r'\.cram$')
_GVCF_PATTERN = re.compile(r'\.g\.vcf(\.gz)?$')
_VCF_PATTERN = re.compile(r'\.vcf(\.gz)?$')

SG_QUERY = graphql.gql(
    """
    a query here to take a project name as a parameter,
    query for all participants in the dataset, and all SGs
    for each SG all analyses of type 'somalier'
    """
)
ANALYSIS_QUERY = graphql.gql(
    """
    a query here to take a project name and list of SGIDs as a parameter,
    query for all analyses of type gVCF, VCF, CRAM in the project for each
    """
)

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
class SomalierDataclass:
    """optional dataclass for Somalier results linked to a sgid."""
    sgid: str
    somalier_file: str | None


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


@cache
def select_somalier_extract_targets(project: str, sgids: list[str]) -> dict[str, str]:
    """
    Take a list of SGIDs, make a parameterised ANALYSIS_QUERY, and process the outputs

    the output of this method should be a dictionary structure:

    {
        sg_id: cram/vcf/gvcf filepath,
        sg_id2: cram/vcf/gvcf filepath,
    }

    > this should re-used the existing logic from JA's somalier implementation, choosing from available data types
    except for one detail - if a short-read SG ID doesn't have a Somalier file, this should throw an error
    maybe the error vs. warning should be behind a config flag, but we expect all short-read CRAMs to have a somalier file

    We also want to remove the possibility of picking up RNA-seq 'cram' analyses
    """
    pass
