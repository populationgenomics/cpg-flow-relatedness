"""
Reading and writing the `somalier_flags` list on a sequencing group's meta in Metamist.

Both the pipeline's reconciler and the resolve CLI write this list, and the mutation replaces it
wholesale rather than patching entries, so a second independently-written copy of it is a second
chance to drop every flag on an SG. It lives here once.
"""

from rd_qc.utils import sg_ids_tag

from metamist.graphql import gql, query

SOMALIER_FLAGS_KEY = 'somalier_flags'

DATASET_SG_META_QUERY = gql(
    """
    query datasetSgMeta($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups {
                id
                meta
            }
        }
    }
    """
)

SG_META_MUTATION = gql(
    """
    mutation updateSgMeta($dataset: String!, $sgId: String!, $sgMeta: JSON!) {
        sequencingGroup {
            updateSequencingGroup(
                project: $dataset
                sequencingGroup: {id: $sgId, meta: $sgMeta}
            ) {
                id
                meta
            }
        }
    }
    """
)


def sequencing_group_key(flag: dict, sg_id: str) -> str:
    """
    Sorted, underscore-joined SG IDs that this flag involves.

    Pairwise flags (self-relatedness, relatedness) are about two SGs but are recorded against
    only the first of the pair, so this key is what lets a reader work out which SGs a flag
    touches without needing per-category knowledge of where the partner ID lives. Per-SG flags
    (sex inference) key on the SG that owns them.
    """
    sg_id_1, sg_id_2 = flag.get('sg_id_1'), flag.get('sg_id_2')
    if sg_id_1 and sg_id_2:
        return sg_ids_tag([sg_id_1, sg_id_2])
    return sg_id


def read_dataset_sg_meta(dataset: str) -> list[dict]:
    """Every sequencing group in the dataset, as `{'id': ..., 'meta': ...}`."""
    response = query(DATASET_SG_META_QUERY, variables={'dataset': dataset})
    return response['project']['sequencingGroups']


def read_sg_flags(dataset: str, sg_id: str) -> list[dict] | None:
    """
    The stored Somalier flags for one sequencing group.

    `None` means the dataset has no such SG, which is a different problem from an SG that has
    never been flagged, and the two deserve different error messages.
    """
    for sg in read_dataset_sg_meta(dataset):
        if sg['id'] == sg_id:
            return (sg.get('meta') or {}).get(SOMALIER_FLAGS_KEY, [])
    return None


def write_sg_flags(dataset: str, sg_id: str, flags: list[dict]) -> None:
    """Replace one sequencing group's whole Somalier flag list."""
    query(
        SG_META_MUTATION,
        variables={'dataset': dataset, 'sgId': sg_id, 'sgMeta': {SOMALIER_FLAGS_KEY: flags}},
    )
