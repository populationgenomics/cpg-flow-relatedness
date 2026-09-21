"""
This script parses "somalier relate" (https://github.com/brentp/somalier) outputs,
and returns a report whether sex and pedigree matches the provided PED file.

Script can send a report to a Slack channel. To enable that, set SLACK_TOKEN
and SLACK_CHANNEL environment variables, and add "Seqr Loader" app into
a channel with:

/invite @Seqr Loader
"""

import contextlib
import json
from argparse import ArgumentParser
from dataclasses import asdict
from typing import Any

import pandas as pd
from loguru import logger
from peddy import Ped, Sample

from rd_qc.utils import (
    DEGREE_IDENTICAL,
    NO_RELATIONSHIP_LABEL,
    UNSPECIFIED_RELATED,
    VERDICT_OK,
    VERDICT_REFINEMENT,
    SomalierFlag,
    SomalierRelatednessFlag,
    SomalierSexInferenceFlag,
    expected_relationship_label,
    get_project_consanguineous_sg_ids,
    infer_degree,
    refine_expected_relationship,
    relatedness_verdict,
)

from cpg_flow.metamist import get_metamist
from cpg_utils import config, slack, to_path
from cpg_utils.metamist_registration import create_new

# Reporting buckets, most serious first
MISMATCH_BUCKETS = ('identical', 'lost_relationship', 'extra_relationship', 'refinement')


_messages: list[str] = []


def info(msg):
    _messages.append(msg)
    logger.info(msg)


def warning(msg):
    _messages.append(msg)
    logger.warning(msg)


def error(msg):
    _messages.append(msg)
    logger.error(msg)


def _format_mismatch_line(s1, s2, expected_ped_s1, expected_ped_s2, expected_rel, inferred_rel, row) -> str:
    fam1 = expected_ped_s1.family_id if expected_ped_s1 else None
    fam2 = expected_ped_s2.family_id if expected_ped_s2 else None
    if fam1 == fam2:
        line = f'{fam1}: {s1} - {s2}'
    else:
        line = s1 + (f' ({fam1})' if fam1 and fam1 != s1 else '')
        line += ' - '
        line += s2 + (f' ({fam2})' if fam2 and fam2 != s2 else '')
    # The unspecified case drops the 'provided:' wrapper, since its label already reads as a
    # statement about the pedigree. Matches how the flags report words the same thing.
    provided = (
        NO_RELATIONSHIP_LABEL
        if expected_rel == UNSPECIFIED_RELATED
        else f'provided: "{expected_relationship_label(expected_rel)}"'
    )
    return (
        f'{line}, '
        f'{provided}, '
        f'inferred: "{inferred_rel}", '
        f'kin={row["relatedness"]}, '
        f'ibs0={row["ibs0"]}, '
        f'ibs2={row["ibs2"]}'
    )


def _check_sex(samples_df: pd.DataFrame) -> dict[str, SomalierSexInferenceFlag]:
    info('*Inferred vs. reported sex:*')
    samples_df.sex = samples_df.sex.apply(lambda x: {1: 'male', 2: 'female'}.get(x, 'unknown'))
    samples_df.original_pedigree_sex = samples_df.original_pedigree_sex.apply(lambda x: {'-9': 'unknown'}.get(x, x))
    samples_df['x_het_ratio'] = samples_df.X_het / samples_df.X_hom_alt.replace(0, 1)  # avoid division by zero
    samples_df['x_depth_ratio'] = 2 * samples_df.X_depth_mean / samples_df.gt_depth_mean
    samples_df['y_depth_ratio'] = 2 * samples_df.Y_depth_mean / samples_df.gt_depth_mean
    missing_inferred_sex = samples_df.sex == 'unknown'
    missing_provided_sex = samples_df.original_pedigree_sex == 'unknown'
    mismatching_female = (samples_df.sex == 'female') & (samples_df.original_pedigree_sex == 'male')
    mismatching_male = (samples_df.sex == 'male') & (samples_df.original_pedigree_sex == 'female')
    mismatching_sex = mismatching_female | mismatching_male
    mismatching_other = (
        (samples_df.sex != samples_df.original_pedigree_sex) & (~mismatching_female) & (~mismatching_male)
    )
    matching_sex = ~mismatching_sex & ~mismatching_other

    sex_mismatches_by_sgid: dict[str, SomalierSexInferenceFlag] = {
        row.sample_id: SomalierSexInferenceFlag(
            category='sex_inference_mismatch',
            provided=str(row.original_pedigree_sex),
            inferred=str(row.sex),
            mean_depth=float(row.gt_depth_mean),
            x_het_ratio=float(row.x_het_ratio),
            x_depth_ratio=float(row.x_depth_ratio),
            y_depth_ratio=float(row.y_depth_ratio),
            x_sites=int(row.X_n),
            p_middling_ab=float(row.p_middling_ab),
        )
        for _, row in samples_df[mismatching_sex].iterrows()
    }

    def _print_stats(df_filter: pd.Series) -> None:
        for _, row_ in samples_df[df_filter].iterrows():
            info(
                f' {row_.sample_id} ('
                f'provided: {row_.original_pedigree_sex}, '
                f'inferred: {row_.sex}, '
                f'X het ratio: {row_.x_het_ratio:.2f})'
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

    return sex_mismatches_by_sgid


CO_PARENTS = 'mom-dad'


def _union_is_recorded_consanguineous(
    expected_rel: str,
    sample_1: Sample | dict,
    sample_2: Sample | dict,
    consanguineous_sgs: set[str],
) -> bool:
    """
    Whether a flagged pair of parents have a recorded consanguineous union in the pedigree.
    """
    if expected_rel != CO_PARENTS or not consanguineous_sgs:
        return False
    kids_1 = {kid.sample_id for kid in getattr(sample_1, 'kids', [])}
    kids_2 = {kid.sample_id for kid in getattr(sample_2, 'kids', [])}
    return bool(kids_1 & kids_2 & consanguineous_sgs)


def _mismatch_bucket(expected_rel: str, measured_rel: str, verdict: str) -> str:
    """Which reporting bucket a flagged pair belongs in, most serious first."""
    if measured_rel == DEGREE_IDENTICAL:
        return 'identical'
    if verdict == VERDICT_REFINEMENT:
        return 'refinement'
    # A conflict is 'lost' when the pedigree expected a closer relationship than was measured,
    # which is the sample-swap shape, and 'extra' when the measurement is closer than expected.
    expected_closer = expected_rel not in ('unrelated', 'mom-dad')
    return 'lost_relationship' if expected_closer else 'extra_relationship'


BUCKET_HEADINGS = {
    'identical': '❗❗ {n} sample pair(s) inferred as identical but recorded as different individuals:',
    'lost_relationship': '❗ {n} sample pair(s) that are recorded as related, but measured as less related:',
    'extra_relationship': '⚠️ {n} sample pair(s) that are recorded as unrelated, but measured as related:',
}


def _report_relatedness_findings(buckets: dict[str, list[str]]) -> None:
    for key, heading in BUCKET_HEADINGS.items():
        pairs = buckets.get(key) or []
        if not pairs:
            continue
        info(heading.format(n=len(pairs)))
        for i, pair in enumerate(pairs):
            info(f' {i + 1}. {pair}')
    if buckets.get('refinement'):
        # Counted but not listed: these are routinely the bulk of the flags and listing them
        # drowns the buckets above, which are the ones that need a decision.
        info(
            f'ℹ️ {len(buckets["refinement"])} pair(s) measured as related where the pedigree records '  # noqa: RUF001
            f'no relationship between them (usually only one parent on file). See the flags report.',
        )
    if not any(buckets.get(key) for key in BUCKET_HEADINGS):
        info('✅ Measured relatedness matches the pedigree for every pair.')
    info('')


def _check_relatedness(
    pairs_df,
    expected_ped: Ped,
    bad_ids: list,
    consanguineous_sgs: set[str] | None = None,
) -> dict[str, list[SomalierRelatednessFlag]]:
    """
    Compare what the pedigree expects against what somalier measured, pair by pair.

    peddy is used only on the expected pedigree. The relationship the data supports
    comes from the kinship coefficient and ibs0 via `infer_degree`.

    `consanguineous_sgs` lets a co-parent pair be excused when the pedigree already records their
    union. Empty or omitted means nothing is excused.
    """
    info('*Relatedness:*')
    expected_ped_sample_by_id: dict[str, Sample] = {s.sample_id: s for s in expected_ped.samples()}

    buckets: dict[str, list[str]] = {key: [] for key in MISMATCH_BUCKETS}

    relatedness_flags_by_sg_id: dict[str, list[SomalierRelatednessFlag]] = {}
    for idx, row in pairs_df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        if s1 in bad_ids or s2 in bad_ids:
            continue

        expected_ped_s1 = expected_ped_sample_by_id.get(s1, {})
        expected_ped_s2 = expected_ped_sample_by_id.get(s2, {})
        with contextlib.redirect_stderr(None), contextlib.redirect_stdout(None):
            if expected_ped_s1 and expected_ped_s2:
                expected_rel = refine_expected_relationship(
                    expected_ped.relation(expected_ped_s1, expected_ped_s2),
                    expected_ped_s1.family_id,
                    expected_ped_s2.family_id,
                )
            else:
                expected_rel = 'unknown'

        measured_rel = infer_degree(row['relatedness'], row['ibs0'], row['n'])
        verdict = relatedness_verdict(expected_rel, measured_rel)

        pairs_df.loc[idx, 'provided_rel'] = expected_rel
        pairs_df.loc[idx, 'inferred_rel'] = measured_rel

        if verdict == VERDICT_OK:
            continue

        if _union_is_recorded_consanguineous(
            expected_rel, expected_ped_s1, expected_ped_s2, consanguineous_sgs or set()
        ):
            logger.info(f'{s1} - {s2}: related co-parents, but a child records the union. Not flagged.')
            continue

        # Make sure that the s1 / s2 sample IDs are sorted to ensure consistent keying
        if s1 > s2:
            s1, s2 = s2, s1
            expected_ped_s1, expected_ped_s2 = expected_ped_s2, expected_ped_s1

        line = _format_mismatch_line(s1, s2, expected_ped_s1, expected_ped_s2, expected_rel, measured_rel, row)
        # peddy .samples() yields Sample objects (attribute access), but the dict lookup above
        # falls back to {} when a sample is missing, so guard both cases with getattr.
        family_external_id = (
            getattr(expected_ped_s1, 'family_id', None) or getattr(expected_ped_s2, 'family_id', None) or 'unknown'
        )
        relatedness_flags_by_sg_id.setdefault(s1, []).append(
            SomalierRelatednessFlag(
                category='relatedness_mismatch',
                sg_id_1=s1,
                sg_id_2=s2,
                family_external_id=family_external_id,
                expected_relationship=expected_rel,
                inferred_relationship=measured_rel,
                relatedness=row['relatedness'],
                ibs0=row['ibs0'],
                ibs2=row['ibs2'],
                verdict=verdict,
            ),
        )

        buckets[_mismatch_bucket(expected_rel, measured_rel, verdict)].append(line)

    _report_relatedness_findings(buckets)

    return relatedness_flags_by_sg_id


def produce_flags(
    somalier_samples: str,
    somalier_pairs: str,
    expected_ped_path: str,
    consanguineous_sgs: set[str] | None = None,
) -> tuple[dict[str, list[SomalierFlag]], pd.DataFrame, pd.DataFrame]:
    """
    Read the somalier relate outputs and produce every flag they imply, keyed by SG id.

    Reads only local files, so callers that have Metamist access pass `consanguineous_sgs` in
    rather than having this query for it. Returns the flags plus both dataframes, which `run`
    needs for its logging.
    """
    logger.info(somalier_samples)
    samples_df = pd.read_csv(somalier_samples, delimiter='\t')
    pairs_df = pd.read_csv(somalier_pairs, delimiter='\t')
    with to_path(expected_ped_path).open() as f:
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

    sex_mismatches_by_sgid = _check_sex(samples_df)
    relatedness_flags_by_sg_id = _check_relatedness(
        pairs_df,
        expected_ped,
        bad_ids,
        consanguineous_sgs,
    )

    all_flags_by_sg_id: dict[str, list[SomalierFlag]] = {}
    for sg_id, flag in sex_mismatches_by_sgid.items():
        all_flags_by_sg_id[sg_id] = [flag]
    for sg_id, flags in relatedness_flags_by_sg_id.items():
        all_flags_by_sg_id.setdefault(sg_id, []).extend(flags)

    return all_flags_by_sg_id, samples_df, pairs_df


def run(
    dataset: str,
    title: str,
    sg_ids: list[str],
    expected_ped: str,
    somalier_pairs: str,
    somalier_samples: str,
    output_pairs: str,
    output_samples: str,
    output_html: str,
    base_output_html: str,
    html_url: str,
    output_json: str,
):
    """Report pedigree inconsistencies given somalier outputs."""

    dataset = get_metamist().get_metamist_proj(dataset)

    all_flags_by_sg_id, samples_df, pairs_df = produce_flags(
        somalier_samples=somalier_samples,
        somalier_pairs=somalier_pairs,
        expected_ped_path=expected_ped,
        consanguineous_sgs=get_project_consanguineous_sg_ids(dataset),
    )

    print_contents(
        samples_df,
        pairs_df,
        somalier_samples,
        somalier_pairs,
    )

    if dataset and html_url:
        title = f'*[{dataset}]* <{html_url}|{title or "Somalier pedigree report"}>'
    elif not title:
        title = 'Somalier pedigree report'
    text = '\n'.join([title, *_messages])

    if config.config_retrieve(['somalier_pedigree', 'send_to_slack'], default=True):
        slack.send_message(text)

    result: dict[str, Any] = {
        'dataset': dataset,
        'html_url': html_url,
        'n_samples_flagged': len([flags for flags in all_flags_by_sg_id.values() if flags]),
        'relatedness_flags': {sg_id: [asdict(flag) for flag in flags] for sg_id, flags in all_flags_by_sg_id.items()},
    }
    if output_json:
        with to_path(output_json).open('w') as f:
            json.dump(result, f, indent=2)

    with to_path(output_samples).open('w') as f:
        f.write(to_path(somalier_samples).read_text())
    with to_path(output_pairs).open('w') as f:
        f.write(to_path(somalier_pairs).read_text())

    # Now create the web analysis for the whole dataset
    create_new(
        project=dataset,
        output=output_html,
        analysis_type='web',
        sgs=sg_ids,
        meta={'stage': 'SomalierPedigreeCheck'},
        secondary={
            'base_html': base_output_html,
            'samples': output_samples,
            'pairs': output_pairs,
            'json': output_json,
        },
    )
    logger.info(f'Registered web analysis for {dataset} at {output_html}')

    return result


def print_contents(
    samples_df,
    pairs_df,
    somalier_samples_fpath,
    somalier_pairs_fpath,
):
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
    parser.add_argument('--dataset', help='Dataset name')
    parser.add_argument('--title', required=True, help='Report title')
    parser.add_argument('--sg-ids', nargs='+', required=True, help='space-separated SG IDs')
    parser.add_argument(
        '--expected-ped',
        required=True,
        help='PED file with expected pedigree',
    )
    parser.add_argument(
        '--somalier-samples',
        required=True,
        help='Somalier samples.tsv file from relate job',
    )
    parser.add_argument(
        '--somalier-pairs',
        required=True,
        help='Somalier pairs.tsv file from relate job',
    )
    parser.add_argument('--output-pairs', required=True, help='gs:// path to output pairs TSV')
    parser.add_argument('--output-samples', required=True, help='gs:// path to output samples TSV')
    parser.add_argument('--output-html', required=True, help='gs:// path to HTML (namespaced by AR GUID)')
    parser.add_argument('--base-output-html', required=True, help='gs:// path to HTML (fixed, not namespaced)')
    parser.add_argument('--html-url', help='Web-accessible HTML path (namespaced by AR GUID)')
    parser.add_argument('--output-json', required=True, help='gs:// path to JSON output for results')
    args = parser.parse_args()
    run(
        dataset=args.dataset,
        title=args.title,
        sg_ids=args.sg_ids,
        expected_ped=args.expected_ped,
        somalier_pairs=args.somalier_pairs,
        somalier_samples=args.somalier_samples,
        output_pairs=args.output_pairs,
        output_samples=args.output_samples,
        output_html=args.output_html,
        base_output_html=args.base_output_html,
        html_url=args.html_url,
        output_json=args.output_json,
    )
