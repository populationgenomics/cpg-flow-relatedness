"""
Render the Somalier flags report from real, locally downloaded check_pedigree outputs.

Reads nothing but local files plus two read-only Metamist queries, and writes nothing anywhere
except the output HTML. In particular it never mutates SG meta and never registers an analysis,
so it is safe to point at a production dataset.

Collect the inputs first. Find the run's AR GUID:

    SM_ENVIRONMENT=production python -c "
    from metamist.graphql import gql, query
    q = gql('query(\\$d: String!){project(name:\\$d){analyses(type:{eq:\\"web\\"},'
            'meta:{stage:\\"SomalierPedigreeCheck\\"}){timestampCompleted outputs}}}')
    print(query(q, variables={'d': 'DATASET'}))"

then pull the two files it points at:

    gcloud storage cp 'gs://cpg-<ds>[-test]/[<subdir>/]somalier_checks/pedigree/<AR_GUID>/*.tsv' \\
        local_data/<ds>/

Only samples.tsv and pairs.tsv are needed, because they are raw `somalier relate` output and so
are unaffected by any change to our own code. The expected PED is rebuilt from Metamist on every
run via build_ped_content, exactly as the stage does, and the flags are re-derived locally from
the TSVs, so the run's checks.json is never consulted unless you ask for it. Then:

    SM_ENVIRONMENT=production uv run python testing_scripts/local_pedigree_report.py \\
        --input-dir local_data/<ds> --dataset <ds> --output local_data/<ds>/report.html

Self-relatedness is out of scope here: this harness only regenerates the sex-inference and
pedigree-relatedness flags that check_pedigree produces. Any self-relatedness flags already in an
SG's meta are passed through untouched rather than reconciled against an empty set, which would
wrongly mark every one of them resolved.
"""

import csv
import json
import sys
from argparse import ArgumentParser
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

REPO_ROOT = Path(__file__).absolute().parent.parent
# src/ on the path because a local .venv often holds a stale non-editable rd_qc (see the
# pythonpath note in pyproject.toml).
sys.path.insert(0, str(REPO_ROOT / 'src'))

from loguru import logger  # noqa: E402

from rd_qc.scripts.check_pedigree import produce_flags  # noqa: E402
from rd_qc.scripts.record_somalier_flags import (  # noqa: E402
    reconcile_sg_somalier_relatedness_flags,
    reconcile_sg_somalier_sex_inference_flags,
    sequencing_group_key,
)
from rd_qc.scripts.somalier_flags_report import (  # noqa: E402
    FLAG_CLASSES,
    SGInfo,
    collect_somalier_flags,
    get_sg_infos,
    group_by_family,
    referenced_sg_ids,
    render_report,
    split_active_resolved,
    summarise_flags,
)
from rd_qc.utils import build_ped_content, get_project_sgs_and_fingerprints  # noqa: E402

from cpg_utils.config import set_config_paths  # noqa: E402
from metamist.graphql import gql, query  # noqa: E402

# PED sex codes, as build_ped_content writes them.
PED_SEX = {'1': 'male', '2': 'female'}

# Deliberately not the report's DATASET_SGS_QUERY: this one takes the sequencing type and
# technology as plain arguments so the harness needs no workflow config to run.
DATASET_SGS_QUERY = gql(
    """
    query datasetSgs($dataset: String!, $seqType: String!, $seqTech: String!) {
        project(name: $dataset) {
            sequencingGroups(type: {eq: $seqType}, technology: {eq: $seqTech}) {
                id
                meta
                type
            }
        }
    }
    """
)


def find_input(input_dir: Path, suffix: str) -> Path:
    """Locate one input by suffix, so the harness does not care what the dataset is called."""
    matches = sorted(p for p in input_dir.iterdir() if p.name.endswith(suffix))
    if not matches:
        raise SystemExit(f'No *{suffix} found in {input_dir}')
    if len(matches) > 1:
        raise SystemExit(f'Ambiguous: {len(matches)} files match *{suffix} in {input_dir}')
    return matches[0]


# ---------------------------------------------------------------------------
# Metamist, read-only, cached to disk so repeat renders need no network
# ---------------------------------------------------------------------------
def load_cache(cache: Path) -> dict:
    """Load the snapshot, preserving any keys this version does not know about."""
    raw = json.loads(cache.read_text()) if cache.exists() else {}
    raw.setdefault('sgs', None)
    raw.setdefault('infos', {})
    raw.setdefault('expected_ped', None)
    return raw


def save_cache(cache: Path, snapshot: dict) -> None:
    cache.write_text(json.dumps(snapshot, indent=2))
    logger.info(f'Updated Metamist snapshot {cache}')


def as_sg_infos(raw: dict) -> dict[str, SGInfo]:
    """Rebuild SGInfo objects from the cached JSON, restoring the fastq pair tuples."""
    return {
        sg_id: SGInfo(**{**fields, 'fastq_pairs': [tuple(pair) for pair in fields['fastq_pairs']]})
        for sg_id, fields in raw.items()
    }


def write_local_config(input_dir: Path, dataset: str, access_level: str, seq_type: str, seq_tech: str) -> Path:
    """
    The minimum CPG config that rd_qc.utils' Metamist helpers need, so this harness can run
    without an analysis-runner invocation. Written into the input directory to stay inspectable.
    """
    path = input_dir / 'local_config.toml'
    path.write_text(
        '[workflow]\n'
        f"dataset = '{dataset}'\n"
        f"access_level = '{access_level}'\n"
        f"sequencing_type = '{seq_type}'\n"
        f"sequencing_technology = '{seq_tech}'\n"
    )
    set_config_paths([str(path)])
    return path


def generate_expected_ped(dataset: str, output: Path) -> str:
    """
    Build the expected PED from Metamist, the same way SomalierPedigreeCheck does at
    orchestration time (stages.py:192), rather than trusting a copy downloaded from GCS.

    filter_sgs=True mirrors the stage, restricting to SGs that match the configured sequencing
    type and technology.
    """
    index = get_project_sgs_and_fingerprints(dataset, filter_sgs=True)
    content = build_ped_content(dataset, index)
    output.write_text(content)
    logger.info(f'Generated expected PED for {len(index.by_sg)} SG(s) from Metamist: {output}')
    return content


def warn_on_stale_provided_sex(ped_content: str, samples_tsv: Path) -> None:
    """
    somalier bakes the PED's sex column into samples.tsv as `original_pedigree_sex`, and that
    frozen value is what `_check_sex` compares against. So regenerating the PED refreshes the
    relatedness half of the check but not the sex half. Warn when the two disagree, which means
    the pedigree changed in Metamist after the relate job ran.
    """
    ped_sex = {}
    for line in ped_content.splitlines():
        parts = line.split('\t')
        if len(parts) >= 5:  # noqa: PLR2004
            ped_sex[parts[1]] = PED_SEX.get(parts[4], 'unknown')

    stale = []
    with samples_tsv.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            provided = row.get('original_pedigree_sex') or 'unknown'
            provided = 'unknown' if provided == '-9' else provided
            fresh = ped_sex.get(row['sample_id'])
            if fresh is not None and fresh != provided:
                stale.append(f'{row["sample_id"]} (somalier saw {provided}, Metamist now says {fresh})')

    if stale:
        logger.warning(
            f'{len(stale)} sample(s) whose provided sex changed in Metamist since the relate job ran. '
            'The relatedness flags below are fresh, but their sex-inference flags are computed against '
            f'the stale value baked into samples.tsv: {"; ".join(stale[:10])}'
        )


def fetch_dataset_sgs(dataset: str, seq_type: str, seq_tech: str) -> list[dict]:
    """Every SG in the dataset matching the sequencing type and technology, with its current meta."""
    response = query(DATASET_SGS_QUERY, variables={'dataset': dataset, 'seqType': seq_type, 'seqTech': seq_tech})
    sgs = response['project']['sequencingGroups']
    logger.info(f'{dataset} :: {len(sgs)} sequencing groups match {seq_type}/{seq_tech}')
    return sgs


# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------
def flags_from_checks_json(path: Path) -> dict[str, list[dict]]:
    """The fast path: reuse the flags check_pedigree already wrote, instead of re-deriving them."""
    payload = json.loads(path.read_text())
    return payload.get('relatedness_flags') or {}


def usable_existing_flags(sg_id: str, current: list[dict]) -> list[dict]:
    """Drop any existing meta flag that no longer fits its dataclass, rather than crashing on it."""
    keep = []
    for flag in current:
        flag_class = FLAG_CLASSES.get((flag or {}).get('category'))
        if flag_class is None:
            logger.warning(f'{sg_id} :: ignoring existing flag with unrecognised category')
            continue
        try:
            flag_class(**flag)
        except TypeError as exc:
            logger.warning(f'{sg_id} :: ignoring unusable existing {flag["category"]} flag: {exc}')
            continue
        keep.append(flag)
    return keep


def reconcile_locally(sg: dict, new_flags: list[dict], today: str) -> list:
    """
    Reconcile fresh flags against whatever is already in the SG's meta, entirely in memory.

    This is reconcile_sg_somalier_flags minus the SG_META_MUTATION, reusing the same three pure
    per-category reconcile functions so the resolved/unresolved lifecycle behaves exactly as it
    would in production.
    """
    sg_id = sg['id']
    current = usable_existing_flags(sg_id, (sg.get('meta') or {}).get('somalier_flags') or [])

    for flag in (*current, *new_flags):
        flag['sequencing_group_key'] = sequencing_group_key(flag, sg_id)

    new_sex = {
        (f['provided'], f['inferred']): f for f in new_flags if f['category'] == 'sex_inference_mismatch'
    }
    new_pedigree = {
        (
            f['sg_id_1'],
            f['sg_id_2'],
            f['family_external_id'],
            f['expected_relationship'],
            f['inferred_relationship'],
        ): f
        for f in new_flags
        if f['category'] == 'relatedness_mismatch'
    }

    final = []
    _, sex_final = reconcile_sg_somalier_sex_inference_flags(sg, new_sex, current, today)
    final.extend(sex_final)
    _, pedigree_final = reconcile_sg_somalier_relatedness_flags(sg, new_pedigree, current, today)
    final.extend(pedigree_final)

    # Self-relatedness is not regenerated here, so carry it through as-is.
    self_class = FLAG_CLASSES['self_relatedness_mismatch']
    final.extend(self_class(**f) for f in current if f['category'] == 'self_relatedness_mismatch')
    return final


def build_sg_meta(dataset_sgs: list[dict], new_flags_by_sg: dict[str, list[dict]], today: str) -> list[dict]:
    """Produce the DATASET_SGS_QUERY-shaped list the report reads, with locally reconciled flags."""
    out = []
    for sg in dataset_sgs:
        final = reconcile_locally(sg, new_flags_by_sg.get(sg['id'], []), today)
        out.append(
            {
                'id': sg['id'],
                'type': sg.get('type'),
                'meta': {'somalier_flags': [asdict(flag) for flag in final]},
            }
        )
    return out


# ---------------------------------------------------------------------------
def parse_args():
    parser = ArgumentParser(description='Render the Somalier flags report from local check_pedigree outputs.')
    parser.add_argument('--input-dir', type=Path, required=True, help='directory holding the downloaded outputs')
    parser.add_argument('--dataset', required=True, help='Metamist project name, e.g. seqr or validation-test')
    parser.add_argument('--output', type=Path, help='HTML output path (default: <input-dir>/report.html)')
    parser.add_argument('--cache', type=Path, help='snapshot JSON (default: <input-dir>/metamist_snapshot.json)')
    parser.add_argument('--sequencing-type', default='genome')
    parser.add_argument('--sequencing-technology', default='short-read')
    parser.add_argument('--access-level', default='full', help="'full'/'standard' use the dataset name as given")
    parser.add_argument(
        '--expected-ped',
        type=Path,
        help='use this PED instead of building one from Metamist (escape hatch for an edited pedigree)',
    )
    parser.add_argument('--refresh', action='store_true', help='re-query Metamist even if the snapshot exists')
    parser.add_argument('--offline', action='store_true', help='never query Metamist; use the cached snapshot only')
    parser.add_argument(
        '--from-checks-json',
        action='store_true',
        help='reuse the flags in <ds>.checks.json instead of re-deriving them from the TSVs',
    )
    return parser.parse_args()


def resolve_expected_ped(input_dir: Path, dataset: str, snapshot: dict, args) -> tuple[Path, str]:
    """
    Get the expected PED, rebuilt from Metamist on every online run so it always reflects the
    current pedigree rather than whatever the relate job happened to run with.

    The content is cached into the snapshot purely so `--offline` re-renders keep working.
    `--expected-ped` overrides with a local file, the escape hatch for testing an edited pedigree.
    """
    if args.expected_ped:
        content = args.expected_ped.read_text()
        logger.info(f'Using the expected PED given on the command line: {args.expected_ped}')
        return args.expected_ped, content

    path = input_dir / f'{dataset}.metamist.expected.ped'
    if args.offline:
        content = snapshot.get('expected_ped')
        if not content:
            raise SystemExit(f'--offline given but no cached expected PED exists for {dataset} yet')
        path.write_text(content)
        logger.info(f'Using the cached Metamist-derived expected PED: {path}')
        return path, content

    content = generate_expected_ped(dataset, path)
    snapshot['expected_ped'] = content
    return path, content


def derive_flags(input_dir: Path, expected_ped: Path, from_checks_json: bool) -> dict[str, list[dict]]:
    """Produce the flags locally from the downloaded somalier outputs."""
    if from_checks_json:
        new_flags_by_sg = flags_from_checks_json(find_input(input_dir, '.checks.json'))
        logger.warning(
            f'Loaded flags for {len(new_flags_by_sg)} SG(s) straight from checks.json. These were '
            'produced by whatever version of check_pedigree ran in that job, so they will not '
            'reflect any local changes to the check. Drop the flag to re-derive them.'
        )
        return new_flags_by_sg

    flags_by_sg, _, _ = produce_flags(
        somalier_samples=str(find_input(input_dir, '.samples.tsv')),
        somalier_pairs=str(find_input(input_dir, '.pairs.tsv')),
        expected_ped_path=str(expected_ped),
    )
    new_flags_by_sg = {sg_id: [asdict(f) for f in flags] for sg_id, flags in flags_by_sg.items()}
    logger.info(f'Derived flags for {len(new_flags_by_sg)} SG(s) from the somalier TSVs')
    return new_flags_by_sg


def resolve_dataset_sgs(cache: Path, snapshot: dict, args) -> list[dict]:
    """The dataset's SGs and their current meta. Read-only, and cached to disk."""
    if snapshot['sgs'] is not None and not args.refresh:
        logger.info(f'Using {len(snapshot["sgs"])} cached sequencing groups from {cache}')
        return snapshot['sgs']
    if args.offline:
        raise SystemExit(f'--offline given but {cache} holds no SG list yet')
    snapshot['sgs'] = fetch_dataset_sgs(args.dataset, args.sequencing_type, args.sequencing_technology)
    save_cache(cache, snapshot)
    return snapshot['sgs']


def main() -> None:
    args = parse_args()
    input_dir: Path = args.input_dir
    output: Path = args.output or input_dir / 'report.html'
    cache: Path = args.cache or input_dir / 'metamist_snapshot.json'
    today = datetime.now(tz=UTC).isoformat(timespec='seconds')

    write_local_config(
        input_dir,
        dataset=args.dataset,
        access_level=args.access_level,
        seq_type=args.sequencing_type,
        seq_tech=args.sequencing_technology,
    )

    snapshot = load_cache(cache)
    expected_ped, ped_content = resolve_expected_ped(input_dir, args.dataset, snapshot, args)
    save_cache(cache, snapshot)

    new_flags_by_sg = derive_flags(input_dir, expected_ped, args.from_checks_json)
    if not args.from_checks_json:
        warn_on_stale_provided_sex(ped_content, find_input(input_dir, '.samples.tsv'))

    resolve_dataset_sgs(cache, snapshot, args)

    # 3. Reconcile in memory. Nothing is written back to Metamist.
    sequencing_groups = build_sg_meta(snapshot['sgs'], new_flags_by_sg, today)
    flagged = [sf for sf in collect_somalier_flags(sequencing_groups) if sf.flags]

    # 4. Only now do we know which SGs need metadata: the reconciled flags reference SGs that the
    #    freshly derived flags do not, because meta can already hold flags from earlier runs.
    referenced = referenced_sg_ids(flagged)
    missing = [sg_id for sg_id in referenced if sg_id not in snapshot['infos']]
    if missing and not args.offline:
        fetched = get_sg_infos(missing)
        snapshot['infos'].update({sg_id: asdict(info) for sg_id, info in fetched.items()})
        save_cache(cache, snapshot)
    infos = as_sg_infos(snapshot['infos'])

    still_missing = [sg_id for sg_id in referenced if sg_id not in infos]
    if still_missing:
        logger.warning(
            f'No SG info for {len(still_missing)} referenced SG(s), so they will group by SG id '
            f'rather than family: {", ".join(still_missing[:10])}'
        )

    # 5. Render, exactly as the real report does.
    groups = group_by_family(flagged, infos)
    active_groups, resolved_groups = split_active_resolved(groups, infos)
    summary = summarise_flags(flagged, total_sgs=len(sequencing_groups), families_affected=len(active_groups))

    html = render_report(args.dataset, active_groups, resolved_groups, summary=summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)

    logger.info(f'summary: {summary}')
    for group in active_groups:
        logger.info(f'  {group.label:<26} {group.count_summary}')
    logger.info(f'Wrote {output} ({len(html):,} bytes)')


if __name__ == '__main__':
    main()
