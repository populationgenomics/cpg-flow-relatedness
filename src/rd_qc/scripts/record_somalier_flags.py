import json
import os
from argparse import ArgumentParser
from dataclasses import asdict
from datetime import UTC, datetime

from loguru import logger

from rd_qc.utils import SomalierRelatednessFlag, SomalierSelfRelatednessFlag, SomalierSexInferenceFlag

from metamist.graphql import gql, query

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


def compare_somalier_sex_inference_flag(current_flag: dict, new_flag: dict) -> bool:
    """
    Returns True if the two flags refer to the same somalier issue and the current flag is still
    unresolved. Flag identity is (provided + inferred).
    """
    return (
        current_flag.get('provided') == new_flag.get('provided')
        and current_flag.get('inferred') == new_flag.get('inferred')
        and not current_flag.get('resolved', False)
    )


def compare_somalier_self_relatedness_flag(current_flag: dict, new_flag: dict) -> bool:
    """
    Returns True if the two flags refer to the same somalier relate issue and the current flag is still
    unresolved. Flag identity is (sg_id_1 + sg_id_2 + participant_external_id + threshold); the measured
    `relatedness`/`ibs0`/`ibs2` are intentionally excluded because they can drift slightly between relate runs.
    """
    return (
        current_flag.get('sg_id_1') == new_flag.get('sg_id_1')
        and current_flag.get('sg_id_2') == new_flag.get('sg_id_2')
        and current_flag.get('participant_external_id') == new_flag.get('participant_external_id')
        and current_flag.get('threshold') == new_flag.get('threshold')
        and not current_flag.get('resolved', False)
    )


def compare_somalier_relatedness_flag(current_flag: dict, new_flag: dict) -> bool:
    """
    Returns True if the two flags refer to the same somalier relate issue and the current flag is still
    unresolved.
    Flag identity is (sg_id_1 + sg_id_2 + family_external_id + expected_relationship + inferred_relationship);
    the measured `relatedness`/`ibs0`/`ibs2` are intentionally excluded because they can drift between relate runs.
    """
    return (
        current_flag.get('sg_id_1') == new_flag.get('sg_id_1')
        and current_flag.get('sg_id_2') == new_flag.get('sg_id_2')
        and current_flag.get('family_external_id') == new_flag.get('family_external_id')
        and current_flag.get('expected_relationship') == new_flag.get('expected_relationship')
        and current_flag.get('inferred_relationship') == new_flag.get('inferred_relationship')
        and not current_flag.get('resolved', False)
    )


def reconcile_sg_somalier_sex_inference_flags(
    sg: dict,
    new_somalier_sex_inference_flags_by_key: dict[tuple[str, str], dict],
    current_somalier_flags: list[dict],
    today: str,
) -> tuple[dict, list[SomalierSexInferenceFlag]]:
    """
    Reconcile current and new Somalier sex inference flags for a single SG and update its meta in Metamist.
    """
    # Keyed by (provided, inferred) for easy lookup
    existing_flags_by_key = {
        (flag['provided'], flag['inferred']): flag
        for flag in current_somalier_flags
        if flag.get('category') == 'sex_inference_mismatch'
    }
    final_flags: list[SomalierSexInferenceFlag] = []
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0}
    sg_id = sg['id']
    report = 'Somalier'
    logger.info(f'{sg_id} :: Found {len(current_somalier_flags)} existing {report} flags. Reconciling.')
    for flag_key, flag in existing_flags_by_key.items():
        if flag_key not in new_somalier_sex_inference_flags_by_key:
            if not flag['resolved']:
                # Previously-flagged issue no longer present: mark resolved
                flag['resolved'] = True
                flag['resolution_date'] = today
                logger.info(f"{sg_id} :: Marking {report} flag '{flag['provided']}-{flag['inferred']}' as resolved.")
                stats['resolved'] += 1
            else:
                # Already resolved and still absent: keep as-is
                logger.debug(f"{sg_id} :: {report} flag '{flag['provided']}-{flag['inferred']}' remains resolved.")
        elif compare_somalier_sex_inference_flag(flag, new_somalier_sex_inference_flags_by_key[flag_key]):
            # Same unresolved issue is still present: refresh the measured value and but keep resolution status.
            # Identity (provided/inferred) is unchanged so this counts as 'retained', not 'updated'.
            new_flag = new_somalier_sex_inference_flags_by_key[flag_key]
            flag |= {
                key: new_flag[key]
                for key in ['mean_depth', 'x_het_ratio', 'x_depth_ratio', 'y_depth_ratio', 'p_middling_ab']
            }
            logger.info(
                f"{sg_id} :: {report} flag '{flag['provided']}-{flag['inferred']}' "
                'remains unresolved (value refreshed).'
            )
            stats['retained'] += 1
        else:
            # Current flag exists in new run but differs (or was resolved and has reappeared):
            # overwrite with new flag data (which sets resolved=False)
            flag.update(new_somalier_sex_inference_flags_by_key[flag_key])
            logger.info(
                f"{sg_id} :: {report} flag '{flag['provided']}-{flag['inferred']}' updated with new information."
            )
            stats['updated'] += 1
        final_flags.append(SomalierSexInferenceFlag(**flag))

    # Add new flags that aren't already represented in current_somalier_flags
    for flag_key, flag in new_somalier_sex_inference_flags_by_key.items():
        if flag_key in existing_flags_by_key:
            continue
        logger.info(f"{sg_id} :: Adding new {report} flag '{flag['provided']}-{flag['inferred']}'.")
        final_flags.append(SomalierSexInferenceFlag(**flag))
        stats['added'] += 1
    return stats, final_flags


def reconcile_sg_somalier_self_relatedness_flags(
    sg: dict,
    new_somalier_self_relatedness_flags_by_key: dict[tuple[str, str, str, float], dict],
    current_somalier_flags: list[dict],
    today: str,
) -> tuple[dict, list[SomalierSelfRelatednessFlag]]:
    """
    Reconcile current and new Somalier self-relatedness flags for a single SG and update its meta in Metamist.
    """
    # Keyed by (sg_id_1, sg_id_2, participant_external_id, threshold) for easy lookup
    existing_flags_by_key = {
        (f['sg_id_1'], f['sg_id_2'], f['participant_external_id'], f['threshold']): f
        for f in current_somalier_flags
        if f.get('category') == 'self_relatedness_mismatch'
    }
    final_flags: list[SomalierSelfRelatednessFlag] = []
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0}
    sg_id = sg['id']
    report = 'Somalier'
    logger.info(f'{sg_id} :: Found {len(current_somalier_flags)} existing {report} flags. Reconciling.')
    for flag_key, flag in existing_flags_by_key.items():
        if flag_key not in new_somalier_self_relatedness_flags_by_key:
            if not flag['resolved']:
                # Previously-flagged issue no longer present: mark resolved
                flag['resolved'] = True
                flag['resolution_date'] = today
                logger.info(f"{sg_id} :: Marking {report} flag '{flag['sg_id_1']}-{flag['sg_id_2']}' as resolved.")
                stats['resolved'] += 1
            else:
                # Already resolved and still absent: keep as-is
                logger.debug(f"{sg_id} :: {report} flag '{flag['sg_id_1']}-{flag['sg_id_2']}' remains resolved.")
        elif compare_somalier_self_relatedness_flag(flag, new_somalier_self_relatedness_flags_by_key[flag_key]):
            # Same unresolved issue is still present: refresh the measured value and
            # but keep resolution status. Identity (sg_id_1/sg_id_2/participant_external_id/threshold)
            # is unchanged so this counts as 'retained', not 'updated'.
            new_flag = new_somalier_self_relatedness_flags_by_key[flag_key]
            flag['relatedness'] = new_flag['relatedness']
            flag['ibs0'] = new_flag['ibs0']
            flag['ibs2'] = new_flag['ibs2']
            logger.info(
                f"{sg_id} :: {report} flag '{flag['sg_id_1']}-{flag['sg_id_2']}' remains unresolved (value refreshed)."
            )
            stats['retained'] += 1
        else:
            # Current flag exists in new run but differs (or was resolved and has reappeared):
            # overwrite with new flag data (which sets resolved=False)
            flag.update(new_somalier_self_relatedness_flags_by_key[flag_key])
            logger.info(f"{sg_id} :: {report} flag '{flag['sg_id_1']}-{flag['sg_id_2']}' updated with new information.")
            stats['updated'] += 1
        final_flags.append(SomalierSelfRelatednessFlag(**flag))

    for flag_key, flag in new_somalier_self_relatedness_flags_by_key.items():
        if flag_key in existing_flags_by_key:
            continue
        logger.info(f"{sg_id} :: Adding new {report} flag '{flag['sg_id_1']}-{flag['sg_id_2']}'.")
        final_flags.append(SomalierSelfRelatednessFlag(**flag))
        stats['added'] += 1

    return stats, final_flags


def reconcile_sg_somalier_relatedness_flags(
    sg: dict,
    new_somalier_relatedness_flags_by_key: dict[tuple[str, str, str, str, str], dict],
    current_somalier_flags: list[dict],
    today: str,
) -> tuple[dict, list[SomalierRelatednessFlag]]:
    """
    Reconcile current and new Somalier relatedness flags for a single SG and update its meta in Metamist.
    """
    # Keyed by (sg_id_1, sg_id_2, family_external_id, expected_relationship, inferred_relationship) for easy lookup
    existing_flags_by_key = {
        (f['sg_id_1'], f['sg_id_2'], f['family_external_id'], f['expected_relationship'], f['inferred_relationship']): f
        for f in current_somalier_flags
        if f.get('category') == 'relatedness_mismatch'
    }
    final_flags: list[SomalierRelatednessFlag] = []
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0}
    sg_id = sg['id']
    report = 'Somalier'
    logger.info(f'{sg_id} :: Found {len(existing_flags_by_key)} existing {report} flags. Reconciling.')
    for flag_key, flag in existing_flags_by_key.items():
        if flag_key not in new_somalier_relatedness_flags_by_key:
            if not flag['resolved']:
                # Previously-flagged issue no longer present: mark resolved
                flag['resolved'] = True
                flag['resolution_date'] = today
                logger.info(f"{sg_id} :: Marking {report} flag '{flag['category']}' as resolved.")
                stats['resolved'] += 1
            else:
                # Already resolved and still absent: keep as-is
                logger.debug(f"{sg_id} :: {report} flag '{flag['category']}' remains resolved.")
        elif compare_somalier_relatedness_flag(flag, new_somalier_relatedness_flags_by_key[flag_key]):
            # Same unresolved issue is still present: refresh the measured value and but keep resolution status.
            # Identity (sg_id_1/sg_id_2/family_external_id/expected_relationship/inferred_relationship)  # noqa: ERA001
            # is unchanged so this counts as 'retained', not 'updated'.
            new_flag = new_somalier_relatedness_flags_by_key[flag_key]
            flag['value'] = new_flag['value']
            logger.info(f"{sg_id} :: {report} flag '{flag['category']}' remains unresolved (value refreshed).")
            stats['retained'] += 1
        else:
            # Current flag exists in new run but differs (or was resolved and has reappeared):
            # overwrite with new flag data (which sets resolved=False)
            flag.update(new_somalier_relatedness_flags_by_key[flag_key])
            logger.info(f"{sg_id} :: {report} flag '{flag['category']}' updated with new information.")
            stats['updated'] += 1
        final_flags.append(SomalierRelatednessFlag(**flag))
    for flag_key, flag in new_somalier_relatedness_flags_by_key.items():
        if flag_key in existing_flags_by_key:
            continue
        logger.info(f"{sg_id} :: Adding new {report} flag '{flag['category']}'.")
        final_flags.append(SomalierRelatednessFlag(**flag))
        stats['added'] += 1

    return stats, final_flags


def reconcile_sg_somalier_flags(
    sg: dict,
    new_flags_by_sg: dict[str, list[dict]],
    dataset: str,
    today: str,
) -> None:
    """
    Reconcile current and new Somalier flags for a single SG and update its meta in Metamist.

    Marks absent flags as resolved, retains unchanged flags, updates differing flags,
    and adds new flags that didn't previously exist.
    """
    sg_id = sg['id']
    report = 'Somalier'
    # Get all the existing QC flags of the specified type for this SG
    somalier_flags_key = 'somalier_flags'

    current_somalier_flags: list[dict] = (sg['meta'] or {}).get(somalier_flags_key, [])
    unresolved_current_flags = [flag for flag in current_somalier_flags if not flag.get('resolved', False)]

    new_somalier_flags: list[dict] = new_flags_by_sg.get(sg_id, [])

    if not current_somalier_flags and not new_somalier_flags:
        logger.info(f'{sg_id} :: No existing or new {report} flags for this SG, skipping.')
        return  # No existing or new Somalier flags for this SG, skip

    if not unresolved_current_flags and not new_somalier_flags:
        logger.info(f'{sg_id} :: No unresolved existing or new {report} flags for this SG, skipping.')
        return  # No unresolved existing or new Somalier flags for this SG, skip

    new_somalier_sex_inference_flags_by_key = {
        (f['provided'], f['inferred']): f for f in new_somalier_flags if f['category'] == 'sex_inference_mismatch'
    }
    new_somalier_self_relatedness_flags_by_key = {
        (f['sg_id_1'], f['sg_id_2'], f['participant_external_id'], f['threshold']): f
        for f in new_somalier_flags
        if f['category'] == 'self_relatedness_mismatch'
    }
    new_somalier_relatedness_flags_by_key = {
        (f['sg_id_1'], f['sg_id_2'], f['family_external_id'], f['expected_relationship'], f['inferred_relationship']): f
        for f in new_somalier_flags
        if f['category'] == 'relatedness_mismatch'
    }

    # Track the final set of flags to be recorded in Metamist, including resolved, retained, updated, and added flags
    final_flags: list[SomalierSexInferenceFlag | SomalierSelfRelatednessFlag | SomalierRelatednessFlag] = []
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0}

    logger.info(f'{sg_id} :: Found {len(current_somalier_flags)} existing {report} flags. Reconciling.')
    if new_somalier_sex_inference_flags_by_key:
        sex_stats, sex_final_flags = reconcile_sg_somalier_sex_inference_flags(
            sg,
            new_somalier_sex_inference_flags_by_key,
            current_somalier_flags,
            today,
        )
        stats = {k: stats[k] + sex_stats.get(k, 0) for k in stats}
        final_flags.extend(sex_final_flags)
    if new_somalier_self_relatedness_flags_by_key:
        self_relatedness_stats, self_relatedness_final_flags = reconcile_sg_somalier_self_relatedness_flags(
            sg,
            new_somalier_self_relatedness_flags_by_key,
            current_somalier_flags,
            today,
        )
        stats = {k: stats[k] + self_relatedness_stats.get(k, 0) for k in stats}
        final_flags.extend(self_relatedness_final_flags)
    if new_somalier_relatedness_flags_by_key:
        relatedness_stats, relatedness_final_flags = reconcile_sg_somalier_relatedness_flags(
            sg,
            new_somalier_relatedness_flags_by_key,
            current_somalier_flags,
            today,
        )
        stats = {k: stats[k] + relatedness_stats.get(k, 0) for k in stats}
        final_flags.extend(relatedness_final_flags)

    # Perform the mutation to update the SG meta
    query(
        SG_META_MUTATION,
        variables={
            'dataset': dataset,
            'sgId': sg_id,
            'sgMeta': {somalier_flags_key: [asdict(flag) for flag in final_flags]},
        },
    )
    logger.info(
        f'{sg_id} :: Recorded {len(final_flags)} {report} flags in Metamist. '
        f'Resolved: {stats["resolved"]}, Retained: {stats["retained"]}, '
        f'Updated: {stats["updated"]}, Added: {stats["added"]}'
    )


def main(
    dataset: str,
    somalier_self_relatedness_json_dir: str,
    somalier_relatedness_json_path: str,
    sg_ids: list[str],
):
    """
    Reads the qc flags JSON file and the SG mapping file, and updates any flagged QC issues in the
    sequencing group meta in Metamist. The mapping file is required to ensure that the SG was in
    scope for the dataset being processed and to avoid updating unrelated SGs.

    If the SG meta already has a 'somalier_flags' key, it will be updated with the new flags, including
    recording resolution information. If the 'somalier_flags' key does not exist, it will be created.
    """
    today = datetime.now(tz=UTC).isoformat(timespec='seconds')

    # Load the Self Relatedness flags from the JSON files
    somalier_self_relatedness_data_by_sg_id = {}
    for filename in os.listdir(somalier_self_relatedness_json_dir):
        if filename.endswith('.json'):
            path = os.path.join(somalier_self_relatedness_json_dir, filename)
            with open(path) as f:
                d = json.load(f)
                somalier_self_relatedness_data_by_sg_id.update(d['self_relatedness_flags'])

    # Load the Relatedness flags from the JSON file
    with open(somalier_relatedness_json_path) as f:
        somalier_relatedness_data = json.load(f)

    # Query the sequencing groups for the given dataset
    response = query(DATASET_SG_META_QUERY, variables={'dataset': dataset})
    sequencing_groups = response['project']['sequencingGroups']

    # Reconcile each sequencing group's QC flags
    new_flags_by_sg: dict[str, list[dict]] = {}
    for sg_id, flags in somalier_self_relatedness_data_by_sg_id.items():
        new_flags_by_sg.setdefault(sg_id, []).extend(flags)
    for sg_id, flags in somalier_relatedness_data['relatedness_flags'].items():
        new_flags_by_sg.setdefault(sg_id, []).extend(flags)

    for sg in sequencing_groups:
        if sg['id'] in sg_ids:
            reconcile_sg_somalier_flags(sg, new_flags_by_sg, dataset, today)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '--dataset',
        required=True,
        help='Dataset name',
    )
    parser.add_argument(
        '--somalier-self-relatedness-json-dir',
        required=True,
        help='Path to directory containing somalier self-relatedness JSON files',
    )
    parser.add_argument(
        '--somalier-relatedness-json-path',
        required=True,
        help='Path to somalier relatedness JSON file',
    )
    parser.add_argument('--sg-ids', nargs='+', required=True, help='space-separated SG IDs')
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        somalier_self_relatedness_json_dir=args.somalier_self_relatedness_json_dir,
        somalier_relatedness_json_path=args.somalier_relatedness_json_path,
        sg_ids=args.sg_ids,
    )
