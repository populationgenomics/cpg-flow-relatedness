"""
Records long-read CRAM QC flags in each sequencing group's meta in Metamist.

Ported from cpg-flow-rdrnaseq's scripts/record_qc_flags.py. Flags go under `{label.lower()}_qc_flags`
('CRAM' -> 'cram_qc_flags', the key short-read CRAM QC uses), reconciled by (section, flag): a flag that is gone
is marked resolved, one still present has its value refreshed, a new one is added. Only SGs in the run's mapping
file are touched.
"""

import json
from argparse import ArgumentParser
from dataclasses import asdict
from datetime import datetime

from loguru import logger

from rd_qc.flag_store import SG_META_MUTATION, read_dataset_sg_meta
from rd_qc.utils import QcFlag

from cpg_utils.config import dataset_for_access_level
from metamist.graphql import query

_QC_FLAG_FIELDS = set(QcFlag.__dataclass_fields__)


def _make_qc_flag(sg_id: str, flag: dict, source: str) -> QcFlag:
    """Construct a QcFlag from a dict, filtering to known fields."""
    known = {k: flag[k] for k in _QC_FLAG_FIELDS if k in flag}
    try:
        return QcFlag(**known)
    except TypeError:
        logger.error(f'{sg_id} :: Cannot construct QcFlag from {source} flag dict: {flag}')
        raise


def compare_qc_flag(current_flag: dict, new_flag: dict) -> bool:
    """True if the two flags are the same QC issue and the current flag is still unresolved."""
    return (
        current_flag.get('flag') == new_flag.get('flag')
        and current_flag.get('comparison') == new_flag.get('comparison')
        and current_flag.get('threshold') == new_flag.get('threshold')
        and current_flag.get('section') == new_flag.get('section')
        and not current_flag.get('resolved', False)
    )


def reconcile_sg_qc_flags(
    sg: dict,
    new_flags_by_sg: dict[str, list[dict]],
    dataset: str,
    label: str,
    today: datetime,
) -> None:
    """Reconcile current and new QC flags for a single SG and update Metamist."""
    sg_id = sg['id']
    report = f'Long-read {label} MultiQC'

    qc_flags_key = f'{label.lower()}_qc_flags'
    current_qc_flags: list[dict] = (sg['meta'] or {}).get(qc_flags_key, [])
    unresolved_current_flags = [flag for flag in current_qc_flags if not flag.get('resolved', False)]
    new_qc_flags: list[dict] = new_flags_by_sg.get(sg_id, [])

    if not unresolved_current_flags and not new_qc_flags:
        logger.info(f'{sg_id} :: No unresolved existing or new {report} flags for this SG, skipping.')
        return

    new_qc_flags_by_key = {(flag['section'], flag['flag']): flag for flag in new_qc_flags}
    existing_flag_keys = {(flag['section'], flag['flag']) for flag in current_qc_flags}

    final_flags: list[QcFlag] = []
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0}

    for flag in current_qc_flags:
        key = (flag['section'], flag['flag'])
        if key not in new_qc_flags_by_key:
            if not flag.get('resolved', False):
                flag['resolved'] = True
                flag['resolution_date'] = today.isoformat(timespec='seconds')
                logger.info(f"{sg_id} :: Marking {report} flag '{flag['flag']}' as resolved.")
                stats['resolved'] += 1
        elif compare_qc_flag(flag, new_qc_flags_by_key[key]):
            new_flag = new_qc_flags_by_key[key]
            flag['value'] = new_flag['value']
            flag['severity'] = new_flag.get('severity', flag.get('severity', 'fail'))
            logger.info(f"{sg_id} :: {report} flag '{flag['flag']}' remains unresolved (value refreshed).")
            stats['retained'] += 1
        else:
            old_severity = flag.get('severity', 'fail')
            flag.update(new_qc_flags_by_key[key])
            transition = f' ({old_severity} -> {flag["severity"]})' if old_severity != flag.get('severity') else ''
            logger.info(f"{sg_id} :: {report} flag '{flag['flag']}' updated{transition}.")
            stats['updated'] += 1
        final_flags.append(_make_qc_flag(sg_id, flag, 'existing'))

    for flag in new_qc_flags:
        if (flag['section'], flag['flag']) in existing_flag_keys:
            continue
        logger.info(f"{sg_id} :: Adding new {report} flag '{flag['flag']}'.")
        final_flags.append(_make_qc_flag(sg_id, flag, 'new'))
        stats['added'] += 1

    query(
        SG_META_MUTATION,
        variables={
            'dataset': dataset,
            'sgId': sg_id,
            'sgMeta': {qc_flags_key: [asdict(flag) for flag in final_flags]},
        },
    )
    logger.info(
        f'{sg_id} :: Recorded {len(final_flags)} {report} flags. '
        f'Resolved: {stats["resolved"]}, Retained: {stats["retained"]}, '
        f'Updated: {stats["updated"]}, Added: {stats["added"]}'
    )


def main(dataset: str, label: str, qc_flags_json_path: str, sg_id_mapping_file: str) -> None:
    """Read the QC flags JSON from check_lrs_qc and update the in-scope SGs' meta in Metamist."""
    today = datetime.now()  # noqa: DTZ005
    project = dataset_for_access_level(dataset)

    with open(sg_id_mapping_file) as f:
        sg_id_map = dict(line.strip().split('\t') for line in f if line.strip())

    with open(qc_flags_json_path) as f:
        new_flags_by_sg = json.load(f)['qc_flags']

    for sg in read_dataset_sg_meta(project):
        if sg['id'] in sg_id_map:
            reconcile_sg_qc_flags(sg, new_flags_by_sg, project, label, today)


if __name__ == '__main__':
    parser = ArgumentParser(description='Record long-read CRAM QC flags in SG meta')
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--label', required=True, help="Flag key label: 'CRAM' -> cram_qc_flags")
    parser.add_argument('--qc-flags-json', required=True, help='Path to the QC flags JSON from check_lrs_qc')
    parser.add_argument('--sequencing-group-ids-map', required=True, help='Path to the SG IDs mapping TSV')
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        label=args.label,
        qc_flags_json_path=args.qc_flags_json,
        sg_id_mapping_file=args.sequencing_group_ids_map,
    )
