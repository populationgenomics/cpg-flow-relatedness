"""
suggested location for any utility methods or constants used across multiple stages
"""

from dataclasses import dataclass
from functools import cache

from metamist import graphql


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

@dataclass
class SomalierDataclass:
    """optional dataclass for Somalier results linked to a sgid."""
    sgid: str
    somalier_file: str | None


@cache
def get_project_sgs_and_fingerprints(project: str):
    """
    Take the project name, make a parameterised SG_QUERY, and process the outputs
    The returned data should preserve:
        - participant ID at the top level
        - each SG ID
        - the path to a somalier file, or None if there isn't one from the SGID

    The output from this method could be a dictionary structure:
    {
        participant_id: {
            sgid: somalier filepath | None,
            sgid: somalier filepath | None,
            sgid: somalier filepath | None
            ...
        },
        participant_id_2: ...
    }
    """
    pass

def find_sgids_without_somalier():
    """
    Takes input from the method above, processes it to find all SGIDs without a somalier extract.
    Should return a set of SGID strings. Could be cached, but the object will be large-ish, so would bloat the cache.
    Unless the cache matches objects on ID instead of contents... in which case it would be super quick.
    """


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
