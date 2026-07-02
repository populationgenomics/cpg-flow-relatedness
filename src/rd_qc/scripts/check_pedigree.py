#!/usr/bin/env python3

"""
This script parses "somalier relate" (https://github.com/brentp/somalier) outputs,
and returns a report whether sex and pedigree matches the provided PED file.

Script can send a report to a Slack channel. To enable that, set SLACK_TOKEN
and SLACK_CHANNEL environment variables, and add "Seqr Loader" app into
a channel with:

/invite @Seqr Loader
"""

import contextlib
from argparse import ArgumentParser

import pandas as pd
from cpg_utils import config, slack, to_path
from loguru import logger
from peddy import Ped

_messages: list[str] = []


def info(msg):
    """
    Record and forward.
    """
    _messages.append(msg)
    logger.info(msg)


def warning(msg):
    """
    Record and forward.
    """
    _messages.append(msg)
    logger.warning(msg)


def error(msg):
    """
    Record and forward.
    """
    _messages.append(msg)
    logger.error(msg)


def run(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    expected_ped_fpath: str,
    title: str,
    html_url: str | None = None,
    dataset: str | None = None,
):
    """Report pedigree inconsistencies, given somalier outputs."""
    logger.info(somalier_samples_fpath)
    samples_df = pd.read_csv(somalier_samples_fpath, delimiter='\t')
    pairs_df = pd.read_csv(somalier_pairs_fpath, delimiter='\t')
    with to_path(somalier_samples_fpath).open() as f:
        inferred_ped = Ped(f)
    with to_path(expected_ped_fpath).open() as f:
        expected_ped = Ped(f)

    bad = samples_df.gt_depth_mean == 0.0
    if bad.any():
        warning(
            f'⚠️ Excluded {len(samples_df[bad])}/{len(samples_df)} samples with zero '
            f'mean GT depth from pedigree/sex checks: {", ".join(samples_df[bad].sample_id)}',
        )
        info('')
    bad_ids = list(samples_df[bad].sample_id)  # for checking in pairs_df
    samples_df = samples_df[~bad]

    info('*Inferred vs. reported sex:*')
    # Rename Ped sex to human-readable tags
    samples_df.sex = samples_df.sex.apply(lambda x: {1: 'male', 2: 'female'}.get(x, 'unknown'))
    samples_df.original_pedigree_sex = samples_df.original_pedigree_sex.apply(lambda x: {'-9': 'unknown'}.get(x, x))
    missing_inferred_sex = samples_df.sex == 'unknown'
    missing_provided_sex = samples_df.original_pedigree_sex == 'unknown'
    mismatching_female = (samples_df.sex == 'female') & (samples_df.original_pedigree_sex == 'male')
    mismatching_male = (samples_df.sex == 'male') & (samples_df.original_pedigree_sex == 'female')
    mismatching_sex = mismatching_female | mismatching_male
    mismatching_other = (
        (samples_df.sex != samples_df.original_pedigree_sex) & (~mismatching_female) & (~mismatching_male)
    )
    matching_sex = ~mismatching_sex & ~mismatching_other

    def _print_stats(df_filter) -> None:
        for _, row_ in samples_df[df_filter].iterrows():
            info(
                f' {row_.sample_id} ('
                f'provided: {row_.original_pedigree_sex}, '
                f'inferred: {row_.sex}, '
                f'mean depth: {row_.gt_depth_mean})',
            )

    if mismatching_sex.any():
        info(f'❗ {len(samples_df[mismatching_sex])}/{len(samples_df)} PED samples with mismatching sex:')
        _print_stats(mismatching_sex)
    if missing_provided_sex.any():
        info(f'⚠️ {len(samples_df[missing_provided_sex])}/{len(samples_df)} samples with missing provided sex:')
        _print_stats(missing_provided_sex)
    if missing_inferred_sex.any():
        info(f'⚠️ {len(samples_df[missing_inferred_sex])}/{len(samples_df)} samples with failed inferred sex:')
        _print_stats(missing_inferred_sex)
    inferred_cnt = len(samples_df[~missing_inferred_sex])
    matching_cnt = len(samples_df[matching_sex])
    info(
        f'✅ Sex inferred for {inferred_cnt}/{len(samples_df)} samples, matching '
        f'for {matching_cnt if matching_cnt != inferred_cnt else "all"} samples.',
    )
    info('')

    info('*Relatedness:*')
    expected_ped_sample_by_id = {s.sample_id: s for s in expected_ped.samples()}
    inferred_ped_sample_by_id = {s.sample_id: s for s in inferred_ped.samples()}

    mismatching_unrelated_to_related = []
    mismatching_related_to_unrelated = []

    for idx, row in pairs_df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        if s1 in bad_ids or s2 in bad_ids:
            continue

        expected_ped_s1 = expected_ped_sample_by_id.get(s1)
        expected_ped_s2 = expected_ped_sample_by_id.get(s2)
        inferred_ped_s1 = inferred_ped_sample_by_id.get(s1)
        inferred_ped_s2 = inferred_ped_sample_by_id.get(s2)
        # Suppressing all logging output from peddy, otherwise it would clutter the logs
        with contextlib.redirect_stderr(None), contextlib.redirect_stdout(None):
            if expected_ped_s1 and expected_ped_s2:
                expected_rel = expected_ped.relation(expected_ped_s1, expected_ped_s2)
            else:
                expected_rel = 'unknown'
            if inferred_ped_s1 and inferred_ped_s2:
                inferred_rel = inferred_ped.relation(inferred_ped_s1, inferred_ped_s2)
            else:
                inferred_rel = 'unknown'

        if inferred_rel != expected_rel:
            # Constructing a line for a report:
            line = ''
            if (fam1 := expected_ped_s1.family_id if expected_ped_s1 else None) == (
                fam2 := expected_ped_s2.family_id if expected_ped_s2 else None
            ):
                line += f'{fam1}: {s1} - {s2}'
            else:
                line += s1 + (f' ({fam1})' if fam1 and fam1 != s1 else '')
                line += ' - '
                line += s2 + (f' ({fam2})' if fam2 and fam2 != s2 else '')
            line = (
                f'{line}, '
                f'provided: "{expected_rel}", '
                f'inferred: "{inferred_rel}", '
                f'kin={row["relatedness"]}, '
                f'ibs0={row["ibs0"]}, '
                f'ibs2={row["ibs2"]}'
            )

            if (expected_rel == 'unknown' and inferred_rel != 'unknown') or (
                expected_rel == 'unrelated' and inferred_rel != 'unrelated'
            ):
                if row['relatedness'] > 0.1:  # noqa: PLR2004
                    mismatching_unrelated_to_related.append(line)
            else:
                mismatching_related_to_unrelated.append(line)

        pairs_df.loc[idx, 'provided_rel'] = expected_rel
        pairs_df.loc[idx, 'inferred_rel'] = inferred_rel

    if mismatching_unrelated_to_related:
        info(
            f'⚠️ Found {len(mismatching_unrelated_to_related)} '
            f'sample pair(s) that are provided as unrelated, are inferred as '
            f'related:',
        )
        for i, pair in enumerate(mismatching_unrelated_to_related):
            info(f' {i + 1}. {pair}')
    if mismatching_related_to_unrelated:
        info(
            f'❗ Found {len(mismatching_related_to_unrelated)} sample pair(s) '
            f'that are provided as related, but inferred as unrelated:',
        )
        for i, pair in enumerate(mismatching_related_to_unrelated):
            info(f' {i + 1}. {pair}')
    if not mismatching_unrelated_to_related and not mismatching_related_to_unrelated:
        info('✅ Inferred pedigree matches for all provided related pairs.')
    info('')

    print_contents(
        samples_df,
        pairs_df,
        somalier_samples_fpath,
        somalier_pairs_fpath,
    )

    # Constructing Slack message
    if dataset and html_url:
        title = f'*[{dataset}]* <{html_url}|{title or "Somalier pedigree report"}>'
    elif not title:
        title = 'Somalier pedigree report'
    text = '\n'.join([title, *_messages])

    if config.config_retrieve(['workflow', 'somalier_pedigree', 'send_to_slack'], default=True):
        slack.send_message(text)

    analysis_api = AnalysisApi()
    for sg_id in sg_ids:
        sg_meta = {
            'check': 'pedigree',
            'sex_match': sg_id not in sex_mismatch_ids,
            'relatedness_issues': [issue for issue in all_issues if sg_id in issue],
        }
        analysis_api.create_analysis(
            project=dataset,
            analysis=Analysis(
                type='somalier_relate',
                status=AnalysisStatus('completed'),
                sequencing_group_ids=[sg_id],
                outputs={
                    'pairs': output_pairs,
                    'samples': output_samples,
                    'html': output_html,
                },
                meta=sg_meta,
            ),
        )
    logger.info(f'Registered somalier_relate analyses for {len(sg_ids)} SGs')


def print_contents(
    samples_df,
    pairs_df,
    somalier_samples_fpath,
    somalier_pairs_fpath,
):
    """
    Print useful information to manually review pedigree check results
    """
    if len(samples_df) < 400:  # noqa: PLR2004
        samples_str = samples_df.to_string()
        logger.info(f'Somalier results, samples (based on {somalier_samples_fpath}):\n{samples_str}\n')
    if len(pairs_df) < 400:  # noqa: PLR2004
        pairs_str = pairs_df[
            [
                '#sample_a',
                'sample_b',
                'relatedness',
                'ibs0',
                'ibs2',
                'n',
                'expected_relatedness',
            ]
        ].to_string()
        logger.info(f'Somalier results, sample pairs (based on {somalier_pairs_fpath}):\n{pairs_str}\n')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '--somalier-samples',
        required=True,
        help='Path to somalier {prefix}.samples.tsv output file',
    )
    parser.add_argument(
        '--somalier-pairs',
        required=True,
        help='Path to somalier {prefix}.pairs.tsv output file',
    )
    parser.add_argument(
        '--ped',
        required=True,
        help='Path to PED file with expected pedigree',
    )
    parser.add_argument('--title', required=True, help='Report title')
    parser.add_argument('--html-url', help='Somalier HTML URL')
    parser.add_argument('--dataset', help='Dataset name')
    args = parser.parse_args()
    run(
        somalier_samples_fpath=args.somalier_samples,
        somalier_pairs_fpath=args.somalier_pairs,
        expected_ped_fpath=args.ped,
        html_url=args.html_url,
        dataset=args.dataset,
        title=args.title,
    )




