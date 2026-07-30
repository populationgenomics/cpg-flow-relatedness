"""
Check somalier self-relatedness results for a single participant.
If any SG pair has relatedness below the threshold, sends a Slack alert.
Registers the relate results in metamist with QC flags.
Always exits 0.
"""

import csv
import json
from argparse import ArgumentParser
from dataclasses import asdict
from typing import Any

from loguru import logger

from rd_qc.utils import SomalierSelfRelatednessFlag

from cpg_flow.metamist import get_metamist
from cpg_utils import config, slack, to_path
from cpg_utils.metamist_registration import create_new


def write_result_json(
    dataset: str,
    participant_external_id: str,
    relatedness_threshold: float,
    html_url: str | None,
    low_relatedness_pairs: list[dict[str, Any]],
    flags_by_sg_id: dict[str, list[SomalierSelfRelatednessFlag]],
    output_json: str | None = None,
):
    """
    Write a JSON file with the self-relatedness check results.
    """
    result: dict[str, Any] = {
        'dataset': dataset,
        'participant_external_id': participant_external_id,
        'relatedness_threshold': relatedness_threshold,
        'html_url': html_url,
        'n_flags': len(low_relatedness_pairs),
        'self_relatedness_flags': {sg_id: [asdict(flag) for flag in flags] for sg_id, flags in flags_by_sg_id.items()},
    }

    if output_json:
        with to_path(output_json).open('w') as f:
            json.dump(result, f, indent=2)


def run(
    pairs_fpath: str,
    participant_external_id: str,
    dataset: str,
    sg_ids: list[str],
    output_pairs: str,
    output_samples: str,
    output_html: str,
    html_url: str,
    output_json: str,
):
    """
    Check somalier self-relatedness results for a single participant.
    If any SG pair has relatedness below the threshold, sends a Slack alert,
    writes the output files, and registers the relate results in Metamist.
    """

    dataset = get_metamist().get_metamist_proj(dataset)
    logger.info(f'Checking self-relatedness for {participant_external_id} in {dataset}')
    relatedness_threshold = config.config_retrieve(
        ['somalier_self_check', 'relatedness_threshold'],
        0.9,
    )
    logger.info(f'Relatedness threshold: {relatedness_threshold}')

    low_relatedness_pairs: list[dict[str, Any]] = []

    try:
        with open(pairs_fpath) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                relatedness = float(row['relatedness'])
                if relatedness < relatedness_threshold:
                    low_relatedness_pairs.append(
                        {
                            'sample_a': row['#sample_a'],
                            'sample_b': row['sample_b'],
                            'relatedness': relatedness,
                            'ibs0': row['ibs0'],
                            'ibs2': row['ibs2'],
                        },
                    )
    except FileNotFoundError:
        logger.warning(f'Pairs file not found: {pairs_fpath} — skipping')
        return

    passed = len(low_relatedness_pairs) == 0
    if passed:
        # All pairs have relatedness above the threshold, exit early
        logger.info(f'{participant_external_id}: All pairs have relatedness >= {relatedness_threshold}')
        write_result_json(
            dataset=dataset,
            participant_external_id=participant_external_id,
            relatedness_threshold=relatedness_threshold,
            html_url=html_url,
            low_relatedness_pairs=[],
            flags_by_sg_id={},
            output_json=output_json,
        )
        return

    flags_by_sg_id: dict[str, list[SomalierSelfRelatednessFlag]] = {}

    header = f'Self-relatedness check failed for participant {participant_external_id}'
    if html_url:
        header = f'<{html_url}|{header}>'
    lines = [
        f'*[{dataset}]* {header}',
        f'Expected relatedness ~1.0 (threshold: {relatedness_threshold}), found:',
    ]
    for pair in low_relatedness_pairs:
        lines.append(
            f'  {pair["sample_a"]} - {pair["sample_b"]}: '
            f'relatedness={pair["relatedness"]}, '
            f'ibs0={pair["ibs0"]}, ibs2={pair["ibs2"]}',
        )

        # Only necessary to register the flag for the first SG in each pair
        if pair['sample_a'] not in flags_by_sg_id:
            flags_by_sg_id[pair['sample_a']] = []
        flags_by_sg_id[pair['sample_a']].append(
            SomalierSelfRelatednessFlag(
                category='self_relatedness_mismatch',
                sg_id_1=pair['sample_a'],
                sg_id_2=pair['sample_b'],
                participant_external_id=participant_external_id,
                threshold=relatedness_threshold,
                relatedness=float(pair['relatedness']),
                ibs0=int(pair['ibs0']),
                ibs2=int(pair['ibs2']),
            )
        )

    text = '\n'.join(lines)
    logger.warning(text)

    if config.config_retrieve(
        ['somalier_self_check', 'send_to_slack'],
        default=True,
    ):
        slack.send_message(text)

    write_result_json(
        dataset=dataset,
        participant_external_id=participant_external_id,
        relatedness_threshold=relatedness_threshold,
        html_url=html_url,
        low_relatedness_pairs=low_relatedness_pairs,
        flags_by_sg_id=flags_by_sg_id,
        output_json=output_json,
    )

    # Register results in metamist
    meta = {
        'stage': 'SomalierSelfCheck',
        'participant_id': participant_external_id,
        'relatedness_threshold': relatedness_threshold,
    }

    create_new(
        project=dataset,
        output=output_pairs,
        analysis_type='somalier_relate',
        sgs=sg_ids,
        meta=meta,
        secondary={'samples': output_samples, 'html': output_html, 'json': output_json},
    )
    logger.info(f'Registered somalier_relate analysis for participant {participant_external_id}')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--pairs-tsv', required=True)
    parser.add_argument('--participant-id', required=True)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--sg-ids', nargs='+', required=True, help='space-separated SG IDs')
    parser.add_argument('--output-pairs', required=True)
    parser.add_argument('--output-samples', required=True)
    parser.add_argument('--output-html', required=True, help='gs:// path to HTML report')
    parser.add_argument('--html-url', required=True, help='Web-accessible URL for HTML report')
    parser.add_argument('--output-json', required=True, help='JSON output file for results')
    args = parser.parse_args()
    run(
        pairs_fpath=args.pairs_tsv,
        participant_external_id=args.participant_id,
        dataset=args.dataset,
        sg_ids=args.sg_ids,
        output_pairs=args.output_pairs,
        output_samples=args.output_samples,
        output_html=args.output_html,
        html_url=args.html_url,
        output_json=args.output_json,
    )
