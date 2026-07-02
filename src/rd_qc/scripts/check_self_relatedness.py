"""
Check somalier self-relatedness results for a single participant.
If any SG pair has kinship below the threshold, sends a Slack alert.
Registers the relate results in metamist with QC flags.
Always exits 0.
"""

import csv
from argparse import ArgumentParser

from cpg_utils import config, slack
from loguru import logger
from metamist.apis import AnalysisApi
from metamist.models import Analysis, AnalysisStatus


def run(
    pairs_fpath: str,
    participant_id: str,
    dataset: str,
    kinship_threshold: float,
    sg_ids: list[str],
    output_pairs: str,
    output_samples: str,
    output_html: str,
):
    logger.info(f'Checking self-relatedness for {participant_id} in {dataset}')
    logger.info(f'Kinship threshold: {kinship_threshold}')

    low_kinship_pairs = []

    try:
        with open(pairs_fpath) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                relatedness = float(row['relatedness'])
                if relatedness < kinship_threshold:
                    low_kinship_pairs.append(
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

    passed = len(low_kinship_pairs) == 0

    if passed:
        logger.info(f'{participant_id}: All pairs have kinship >= {kinship_threshold}')
    else:
        lines = [
            f'*[{dataset}] Self-relatedness check failed for participant {participant_id}*',
            f'Expected kinship ~1.0 (threshold: {kinship_threshold}), found:',
        ]
        for pair in low_kinship_pairs:
            lines.append(
                f'  {pair["sample_a"]} - {pair["sample_b"]}: '
                f'kinship={pair["relatedness"]}, '
                f'ibs0={pair["ibs0"]}, ibs2={pair["ibs2"]}',
            )

        text = '\n'.join(lines)
        logger.warning(text)

        if config.config_retrieve(
            ['workflow', 'somalier_self_check', 'send_to_slack'],
            default=True,
        ):
            slack.send_message(text)

    # Register results in metamist
    meta = {
        'check': 'identity',
        'participant_id': participant_id,
        'passed': passed,
        'kinship_threshold': kinship_threshold,
        'flagged_pairs': low_kinship_pairs,
    }

    AnalysisApi().create_analysis(
        project=dataset,
        analysis=Analysis(
            type='somalier_relate',
            status=AnalysisStatus('completed'),
            sequencing_group_ids=sg_ids,
            outputs={
                'pairs': output_pairs,
                'samples': output_samples,
                'html': output_html,
            },
            meta=meta,
        ),
    )
    logger.info(f'Registered somalier_relate analysis for participant {participant_id}')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--pairs-tsv', required=True)
    parser.add_argument('--participant-id', required=True)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--kinship-threshold', type=float, default=0.9)
    parser.add_argument('--sg-ids', required=True, help='Comma-separated SG IDs')
    parser.add_argument('--output-pairs', required=True)
    parser.add_argument('--output-samples', required=True)
    parser.add_argument('--output-html', required=True)
    args = parser.parse_args()
    run(
        pairs_fpath=args.pairs_tsv,
        participant_id=args.participant_id,
        dataset=args.dataset,
        kinship_threshold=args.kinship_threshold,
        sg_ids=args.sg_ids.split(','),
        output_pairs=args.output_pairs,
        output_samples=args.output_samples,
        output_html=args.output_html,
    )