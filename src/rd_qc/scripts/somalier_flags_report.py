"""
Queries Metamist for all Somalier flags across a dataset's sequencing groups
and renders them into a report using the somalier_flags_overview.html.jinja template.

The page groups flags by family, because a family is the natural review unit for relatedness and
the only grouping under which a pedigree mismatch and a sex mismatch on the same individual appear
side by side. A pedigree pair that straddles two families is rendered under both of them (the
classic cross-family sample swap, which both families' reviewers need to see) but counted once,
via the per-flag identity key.

See docs/superpowers/specs/2026-09-14-somalier-relatedness-report-design.md.
"""

from argparse import ArgumentParser
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter

import jinja2
from loguru import logger

from rd_qc.utils import (
    SomalierFlag,
    SomalierRelatednessFlag,
    SomalierSelfRelatednessFlag,
    SomalierSexInferenceFlag,
    sg_ids_tag,
)

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, dataset_for_access_level
from cpg_utils.metamist_registration import create_new
from metamist.graphql import gql, query

STAGE_NAME = 'GenerateSomalierFlagsReport'
JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'

# Written as escapes rather than literals so that ruff's RUF001 ambiguous-unicode check stays quiet.
PAIR_SEP = ' ↔ '  # left-right arrow, joins the two members of a pair
DOT_SEP = ' · '  # middle dot, separates identifier parts
CROSS_FAMILY_MARK = '↗'  # north-east arrow, marks a cross-family row
DASH = '—'  # em dash, stands in for a missing value

FLAG_CLASSES: dict[str, type[SomalierFlag]] = {
    'sex_inference_mismatch': SomalierSexInferenceFlag,
    'self_relatedness_mismatch': SomalierSelfRelatednessFlag,
    'relatedness_mismatch': SomalierRelatednessFlag,
}

# Short keys drive the filter chips and the template's data-* attributes.
CATEGORY_KEYS = {
    'sex_inference_mismatch': 'sex',
    'self_relatedness_mismatch': 'self',
    'relatedness_mismatch': 'pedigree',
}

# Hand-maintained because the raw category strings are not collaborator-friendly and no friendly
# name exists anywhere upstream. Small and stable, so a dict is enough, with a fallback to the key.
CATEGORY_LABELS = {
    'sex': 'Sex inference',
    'self': 'Self-relatedness',
    'pedigree': 'Pedigree relatedness',
}

CATEGORY_ORDER = ['sex', 'self', 'pedigree']

# Only render the filter bar when there is actually something to filter.
MIN_GROUPS_FOR_FILTER_BAR = 5

# How many flag lines a family group shows inline before collapsing the rest behind its expand.
# Real datasets have families with 50+ pedigree mismatches (perth-neuro family 8718 has 57), and
# rendering those inline buries every other family under a wall of text.
INLINE_FLAG_LIMIT = 5

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

SGS_INFO_QUERY = gql(
    """
    query sgInfo($sgIds: [String!]!) {
        sequencingGroups(id: {in_: $sgIds}) {
            id
            meta
            type
            technology
            platform
            assays {
                id
                meta
            }
            sample {
                id
                externalIds
                type
                participant {
                    id
                    externalIds
                    families {
                        id
                        externalIds
                    }
                }
            }
        }
    }
    """
)

EXISTING_ANALYSES_QUERY = gql(
    """
    query existingAnalyses($dataset: String!, $metaFilter: JSON!) {
        project(name: $dataset) {
            analyses(type: {eq: "web"}, meta: $metaFilter) {
                id
                outputs
                timestampCompleted
                meta
            }
        }
    }
    """
)


@dataclass
class SGInfo:
    sg_id: str
    sg_type: str
    sg_technology: str
    sg_platform: str
    crams: list[str]
    fastq_pairs: list[tuple[str, str]]
    other_reads: list[str]
    sample_external_id: str
    sample_type: str
    participant_external_id: str
    family_external_id: str


@dataclass(frozen=True)
class SgFlags:
    """The Somalier flags recorded against one sequencing group's meta."""

    sg_id: str
    flags: tuple[SomalierFlag, ...]


@dataclass(frozen=True)
class FlagRow:
    """One display-ready flag line. Strings are built here, since Jinja autoescape is on."""

    category: str
    category_key: str
    category_label: str
    identity: tuple
    sg_key: str
    subject: str
    subject_detail: str
    result: str
    details: tuple[tuple[str, str], ...]
    cross_family: str | None
    resolved: bool
    date_short: str
    date_full: str
    resolution_date_short: str
    resolution_date_full: str
    search_blob: str


@dataclass(frozen=True)
class FamilyGroup:
    """A family (or participant/SG fallback) and every flag row that belongs to it."""

    key: str
    label: str
    flags: tuple[FlagRow, ...]
    sg_infos: tuple[SGInfo, ...]
    counts: dict[str, int]
    count_summary: str
    search_blob: str

    @property
    def total(self) -> int:
        return sum(self.counts.values())


def _group_sort_key(group: FamilyGroup) -> tuple[int, str]:
    """Worst families first, then alphabetical."""
    return (-group.total, group.label)


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------
def _fmt_num(value: object) -> str:
    """
    Raw floats out of Somalier are ugly, so: integers stay integers, values at or above 1 get two
    decimal places, values below 1 get two significant figures, trailing zeros are stripped.

    0.616722 -> '0.62', 64.579124 -> '64.58', 4821.0 -> '4821'.
    """
    if value is None or value == '':
        return DASH
    if isinstance(value, (bool, str)):
        return str(value)
    try:
        num = float(value)
    except (TypeError, ValueError):
        return str(value)
    if num.is_integer():
        return str(int(num))
    text = f'{num:.2f}' if abs(num) >= 1 else f'{num:.2g}'
    return text.rstrip('0').rstrip('.') if '.' in text else text


def _date_parts(value: str | None) -> tuple[str, str]:
    """Return (YYYY-MM-DD for display, full ISO timestamp for the hover title)."""
    if not value:
        return '', ''
    return value[:10], value


def _primary_external_id(obj: dict | None) -> str:
    """Metamist keys the primary external ID under ''. Fall back gracefully."""
    ext = (obj or {}).get('externalIds') or {}
    return ext.get('') or next(iter(ext.values()), '')


def _basename(entry: dict | str | None) -> str | None:
    """Pull a file basename from a reads entry (dict or path string)."""
    if isinstance(entry, dict):
        return entry.get('basename') or ((entry.get('location') or '').rsplit('/', 1)[-1] or None)
    if isinstance(entry, str):
        return entry.rsplit('/', 1)[-1] or None
    return None


def _extract_reads(assays: list[dict]) -> tuple[list[str], list[tuple[str, str]], list[str]]:
    """Group an SG's assay read files into (crams, fastq_pairs, other).

    Each assay's ``meta.reads`` is a list of file entries; ``meta.reads_type``
    tells us whether they're fastq (R1/R2 pairs) or aligned (bam/cram).
    """
    crams: list[str] = []
    fastq_pairs: list[tuple[str, str]] = []
    other: list[str] = []

    for assay in assays:
        meta = assay.get('meta') or {}
        reads = meta.get('reads')
        reads_type = (meta.get('reads_type') or '').lower()
        entries = reads if isinstance(reads, list) else ([reads] if reads else [])
        names = [n for n in (_basename(e) for e in entries) if n]
        if not names:
            continue

        if reads_type == 'fastq':
            # A fastq assay is one (or more) R1/R2 pair(s); pair sequentially.
            for i in range(0, len(names), 2):
                r2 = names[i + 1] if i + 1 < len(names) else ''
                fastq_pairs.append((names[i], r2))
        elif reads_type in ('bam', 'cram'):
            crams.extend(names)
        else:
            # Unknown reads_type, so classify by extension.
            for name in names:
                low = name.lower()
                if low.endswith(('.cram', '.bam')):
                    crams.append(name)
                elif low.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')):
                    fastq_pairs.append((name, ''))
                else:
                    other.append(name)

    return crams, fastq_pairs, other


# ---------------------------------------------------------------------------
# Metamist querying
# ---------------------------------------------------------------------------
def get_sg_infos(sg_ids: list[str]) -> dict[str, SGInfo]:
    """Query Metamist for detailed SG info, keyed by SG id."""
    if not sg_ids:
        logger.warning('No SGs flagged')
        return {}
    logger.info(f'Querying Metamist for detailed info on {len(sg_ids)} SG(s): {", ".join(sg_ids)}')
    started = perf_counter()
    response = query(SGS_INFO_QUERY, variables={'sgIds': sg_ids})
    logger.info(f'Received SG info for {len(sg_ids)} SG(s) in {perf_counter() - started:.1f}s')
    infos: dict[str, SGInfo] = {}
    for sg in response['sequencingGroups']:
        sample = sg.get('sample') or {}
        participant = sample.get('participant') or {}
        families = participant.get('families') or []
        family = families[0] if families else None

        crams, fastq_pairs, other = _extract_reads(sg.get('assays', []))

        infos[sg['id']] = SGInfo(
            sg_id=sg['id'],
            sg_type=sg.get('type') or '',
            sg_technology=sg.get('technology') or '',
            sg_platform=sg.get('platform') or '',
            crams=crams,
            fastq_pairs=fastq_pairs,
            other_reads=other,
            sample_external_id=_primary_external_id(sample),
            sample_type=sample.get('type') or '',
            participant_external_id=_primary_external_id(participant),
            family_external_id=_primary_external_id(family),
        )
    return infos


def get_previous_analysis(dataset: str, meta_filter: dict) -> dict | None:
    """Query Metamist for existing web analyses matching a meta filter and take the most recent one, if any."""
    response = query(EXISTING_ANALYSES_QUERY, variables={'dataset': dataset, 'metaFilter': meta_filter})
    existing_analyses = response['project']['analyses']
    if not existing_analyses:
        return None
    existing_analyses.sort(key=lambda a: a.get('timestampCompleted') or '', reverse=True)
    if len(existing_analyses) == 1:
        # Only one analysis is the current one, so there are no previous analyses to compare against.
        return None
    previous_analysis = existing_analyses[1]  # The second most recent analysis is the previous one
    logger.info(f'Found previous analysis {previous_analysis["id"]} from {previous_analysis["timestampCompleted"]}')
    if not previous_analysis['meta'].get('summary'):
        logger.warning(f'Previous analysis {previous_analysis["id"]} has no summary in meta; skipping')
        return None
    return previous_analysis


# ---------------------------------------------------------------------------
# Reading flags out of meta
# ---------------------------------------------------------------------------
def collect_somalier_flags(sequencing_groups: list[dict]) -> list[SgFlags]:
    """
    Instantiate each SG's Somalier flags, dispatching on ``category`` to the right dataclass.

    A malformed or unrecognised flag is logged and skipped rather than raised: a report that dies
    on one bad meta entry is worse than one that renders the other ninety-nine.
    """
    collected: list[SgFlags] = []
    for sg in sequencing_groups:
        meta = sg.get('meta') or {}
        flags: list[SomalierFlag] = []
        for raw in meta.get('somalier_flags') or []:
            category = (raw or {}).get('category')
            flag_class = FLAG_CLASSES.get(category)
            if flag_class is None:
                logger.warning(f'{sg["id"]} :: skipping Somalier flag with unrecognised category {category!r}')
                continue
            try:
                flags.append(flag_class(**raw))
            except TypeError as exc:
                logger.warning(f'{sg["id"]} :: skipping malformed {category} flag: {exc}')
        collected.append(SgFlags(sg_id=sg['id'], flags=tuple(flags)))
    return collected


def flag_sg_key(flag: SomalierFlag, owning_sg_id: str) -> str:
    """
    The SG key recorded on the flag, derived on the fly when absent.

    ``sequencing_group_key`` is written by record_somalier_flags.py, but flags recorded before that
    field existed do not carry it, so fall back to the category-specific ID fields. This is the
    only place the report reads sg_id_1/sg_id_2 directly.
    """
    if flag.sequencing_group_key:
        return flag.sequencing_group_key
    sg_id_1 = getattr(flag, 'sg_id_1', None)
    sg_id_2 = getattr(flag, 'sg_id_2', None)
    if sg_id_1 and sg_id_2:
        return sg_ids_tag([sg_id_1, sg_id_2])
    return owning_sg_id


def flag_identity(flag: SomalierFlag, sg_key: str) -> tuple:
    """
    A flag's dataset-wide identity, excluding the measured values.

    Mirrors the identity keys reconciliation uses (record_somalier_flags.py), so the report and the
    resolved/unresolved lifecycle agree on what counts as the same flag. Measured values drift
    between relate runs, so they are excluded here as they are there.
    """
    category_key = CATEGORY_KEYS.get(flag.category, flag.category or '')
    if category_key == 'sex':
        return (category_key, sg_key, flag.provided, flag.inferred)
    if category_key == 'self':
        return (category_key, sg_key, flag.participant_external_id, flag.threshold)
    if category_key == 'pedigree':
        return (category_key, sg_key, flag.expected_relationship, flag.inferred_relationship)
    return (category_key, sg_key)


def referenced_sg_ids(sg_flags: list[SgFlags]) -> list[str]:
    """
    Every SG ID any flag touches, including the far member of a pairwise flag.

    Pairwise flags are recorded against the first SG of the pair only, so without this the
    partner's family and participant IDs are never fetched and the pair cannot be grouped.
    """
    sg_ids: set[str] = set()
    for sf in sg_flags:
        if sf.flags:
            sg_ids.add(sf.sg_id)
        for flag in sf.flags:
            sg_ids.update(part for part in flag_sg_key(flag, sf.sg_id).split('_') if part)
    return sorted(sg_ids)


# ---------------------------------------------------------------------------
# Grouping by family
# ---------------------------------------------------------------------------
def _family_group(sg_id: str, infos: dict[str, SGInfo], fallback_participant: str = '') -> tuple[str, str]:
    """Resolve one SG to its (group key, display label), falling back to participant then SG."""
    info = infos.get(sg_id)
    family = (info.family_external_id if info else '') or ''
    if family:
        return family, family
    participant = (info.participant_external_id if info else '') or fallback_participant
    if participant:
        return f'participant:{participant}', f'(no family){DOT_SEP}{participant}'
    return f'sg:{sg_id}', f'(no family){DOT_SEP}{sg_id}'


def _pedigree_family_group(sg_id: str, infos: dict[str, SGInfo], recorded_family: str) -> tuple[str, str]:
    """
    As _family_group, but falling back to the family the flag itself recorded.

    Metamist is preferred because the flag's own ``family_external_id`` is lossy: check_pedigree.py
    collapses it to ``fam1 or fam2 or 'unknown'``, discarding the fact that a pair straddled two
    families at all.
    """
    key, label = _family_group(sg_id, infos)
    if key.startswith(('participant:', 'sg:')) and recorded_family and recorded_family != 'unknown':
        return recorded_family, recorded_family
    return key, label


def _participant_of(sg_id: str, infos: dict[str, SGInfo]) -> str:
    info = infos.get(sg_id)
    return (info.participant_external_id if info else '') or ''


def _group_targets(flag: SomalierFlag, owning_sg_id: str, infos: dict[str, SGInfo]) -> list[tuple[str, str]]:
    """
    The one or two family groups this flag belongs in.

    Two only for a pedigree pair whose members resolve to different families, which is the
    cross-family sample-swap case.
    """
    category_key = CATEGORY_KEYS.get(flag.category)
    if category_key == 'pedigree':
        recorded = flag.family_external_id or ''
        first = _pedigree_family_group(flag.sg_id_1, infos, recorded)
        second = _pedigree_family_group(flag.sg_id_2, infos, recorded)
        return [first] if first[0] == second[0] else [first, second]
    if category_key == 'self':
        return [_family_group(flag.sg_id_1, infos, fallback_participant=flag.participant_external_id)]
    return [_family_group(owning_sg_id, infos)]


def _sex_row_parts(flag: SomalierSexInferenceFlag, owning_sg_id: str, infos: dict[str, SGInfo]) -> dict:
    participant = _participant_of(owning_sg_id, infos)
    return {
        'subject': f'{participant}{DOT_SEP}{owning_sg_id}' if participant else owning_sg_id,
        'subject_detail': '',
        'result': f'provided {flag.provided} / inferred {flag.inferred}',
        'details': (
            ('Mean depth', _fmt_num(flag.mean_depth)),
            ('X het ratio', _fmt_num(flag.x_het_ratio)),
            ('X depth ratio', _fmt_num(flag.x_depth_ratio)),
            ('Y depth ratio', _fmt_num(flag.y_depth_ratio)),
            ('X sites', _fmt_num(flag.x_sites)),
            ('p middling AB', _fmt_num(flag.p_middling_ab)),
        ),
        'search_extra': f'{participant} {owning_sg_id} {flag.provided} {flag.inferred}',
    }


def _self_row_parts(flag: SomalierSelfRelatednessFlag, infos: dict[str, SGInfo]) -> dict:
    return {
        'subject': f'{flag.sg_id_1}{PAIR_SEP}{flag.sg_id_2}',
        'subject_detail': flag.participant_external_id,
        'result': (
            f'relatedness {_fmt_num(flag.relatedness)} (expected ~1.0, threshold {_fmt_num(flag.threshold)})'
        ),
        'details': (
            ('Relatedness', _fmt_num(flag.relatedness)),
            ('Threshold', _fmt_num(flag.threshold)),
            ('IBS0', _fmt_num(flag.ibs0)),
            ('IBS2', _fmt_num(flag.ibs2)),
            ('Participant', flag.participant_external_id or DASH),
        ),
        'search_extra': ' '.join(
            [
                flag.sg_id_1,
                flag.sg_id_2,
                flag.participant_external_id or '',
                _participant_of(flag.sg_id_1, infos),
                _participant_of(flag.sg_id_2, infos),
            ]
        ),
    }


def _pedigree_row_parts(flag: SomalierRelatednessFlag, infos: dict[str, SGInfo]) -> dict:
    participant_1 = _participant_of(flag.sg_id_1, infos)
    participant_2 = _participant_of(flag.sg_id_2, infos)
    detail = f'{participant_1}{PAIR_SEP}{participant_2}' if participant_1 and participant_2 else ''
    return {
        'subject': f'{flag.sg_id_1}{PAIR_SEP}{flag.sg_id_2}',
        'subject_detail': detail,
        'result': f'expected {flag.expected_relationship} / inferred {flag.inferred_relationship}',
        'details': (
            ('Relatedness', _fmt_num(flag.relatedness)),
            ('IBS0', _fmt_num(flag.ibs0)),
            ('IBS2', _fmt_num(flag.ibs2)),
            ('Expected', flag.expected_relationship or DASH),
            ('Inferred', flag.inferred_relationship or DASH),
            ('Family (recorded)', flag.family_external_id or DASH),
        ),
        'search_extra': ' '.join(
            [
                flag.sg_id_1,
                flag.sg_id_2,
                participant_1,
                participant_2,
                flag.family_external_id or '',
                flag.expected_relationship or '',
                flag.inferred_relationship or '',
            ]
        ),
    }


def _flag_to_row(
    flag: SomalierFlag,
    owning_sg_id: str,
    infos: dict[str, SGInfo],
    cross_family: str | None = None,
) -> FlagRow:
    """Flatten one flag into a display-ready row, branching on category for labels and result text."""
    category_key = CATEGORY_KEYS.get(flag.category, flag.category or 'unknown')
    if category_key == 'sex':
        parts = _sex_row_parts(flag, owning_sg_id, infos)
    elif category_key == 'self':
        parts = _self_row_parts(flag, infos)
    else:
        parts = _pedigree_row_parts(flag, infos)

    date_short, date_full = _date_parts(flag.date)
    resolution_short, resolution_full = _date_parts(flag.resolution_date)
    sg_key = flag_sg_key(flag, owning_sg_id)

    return FlagRow(
        category=flag.category or '',
        category_key=category_key,
        category_label=CATEGORY_LABELS.get(category_key, category_key),
        identity=flag_identity(flag, sg_key),
        sg_key=sg_key,
        subject=parts['subject'],
        subject_detail=parts['subject_detail'],
        result=parts['result'],
        details=parts['details'],
        cross_family=cross_family,
        resolved=flag.resolved,
        date_short=date_short,
        date_full=date_full,
        resolution_date_short=resolution_short,
        resolution_date_full=resolution_full,
        search_blob=parts['search_extra'].lower(),
    )


def _count_summary(counts: dict[str, int]) -> str:
    """'3 flags (1 sex inference, 2 pedigree relatedness)'."""
    total = sum(counts.values())
    parts = [f'{counts[key]} {CATEGORY_LABELS.get(key, key).lower()}' for key in CATEGORY_ORDER if counts.get(key)]
    label = 'flag' if total == 1 else 'flags'
    return f'{total} {label} ({", ".join(parts)})' if parts else f'{total} {label}'


def _row_sort_key(row: FlagRow) -> tuple[int, str]:
    order = CATEGORY_ORDER.index(row.category_key) if row.category_key in CATEGORY_ORDER else len(CATEGORY_ORDER)
    return (order, row.subject)


def _build_group(key: str, label: str, rows: list[FlagRow], infos: dict[str, SGInfo]) -> FamilyGroup:
    """Assemble a FamilyGroup, recomputing counts, summary and search text from the given rows."""
    counts: dict[str, int] = {}
    for row in rows:
        counts[row.category_key] = counts.get(row.category_key, 0) + 1

    ordered = sorted(rows, key=_row_sort_key)

    sg_ids: list[str] = []
    for row in ordered:
        for part in row.sg_key.split('_'):
            if part and part not in sg_ids:
                sg_ids.append(part)
    group_infos = tuple(infos[sg_id] for sg_id in sg_ids if sg_id in infos)

    search_parts = [label.lower(), *(row.search_blob for row in ordered)]
    for info in group_infos:
        search_parts.extend(
            [
                info.sample_external_id.lower(),
                info.participant_external_id.lower(),
                info.family_external_id.lower(),
                info.sg_id.lower(),
            ]
        )

    return FamilyGroup(
        key=key,
        label=label,
        flags=tuple(ordered),
        sg_infos=group_infos,
        counts=counts,
        count_summary=_count_summary(counts),
        search_blob=' '.join(part for part in search_parts if part),
    )


def group_by_family(sg_flags: list[SgFlags], infos: dict[str, SGInfo]) -> list[FamilyGroup]:
    """
    Bucket every flag into one family group, or two for a cross-family pedigree pair.

    The duplication for cross-family pairs is deliberate: both families' reviewers need to see a
    suspected swap between them. ``FlagRow.identity`` is what stops it being counted twice.
    """
    buckets: dict[str, tuple[str, list[FlagRow]]] = {}
    for sf in sg_flags:
        for flag in sf.flags:
            targets = _group_targets(flag, sf.sg_id, infos)
            for index, (key, label) in enumerate(targets):
                other = targets[1 - index][1] if len(targets) > 1 else None
                row = _flag_to_row(flag, sf.sg_id, infos, cross_family=other)
                buckets.setdefault(key, (label, []))[1].append(row)

    groups = [_build_group(key, label, rows, infos) for key, (label, rows) in buckets.items()]
    return sorted(groups, key=_group_sort_key)


def split_active_resolved(
    groups: list[FamilyGroup],
    infos: dict[str, SGInfo],
) -> tuple[list[FamilyGroup], list[FamilyGroup]]:
    """
    Rebuild the groups twice over, once holding only unresolved flags and once only resolved ones.

    Counts and summaries are recomputed against the filtered flag set, so a family with both kinds
    appears in both sections showing only the flags relevant to each. A side with nothing left is
    dropped.
    """
    active: list[FamilyGroup] = []
    resolved: list[FamilyGroup] = []
    for group in groups:
        for bucket, wanted in ((active, False), (resolved, True)):
            rows = [row for row in group.flags if row.resolved is wanted]
            if rows:
                bucket.append(_build_group(group.key, group.label, rows, infos))
    return sorted(active, key=_group_sort_key), sorted(resolved, key=_group_sort_key)


# ---------------------------------------------------------------------------
# Summary and filter chips
# ---------------------------------------------------------------------------
def summarise_flags(sg_flags: list[SgFlags], total_sgs: int, families_affected: int) -> dict:
    """
    Dataset-wide, flag-centric counts for the header cards.

    Deduplicates on flag identity before counting. Flags are currently recorded against one SG of a
    pair only, so nothing is duplicated today, but this keeps the counts right if that ever changes.
    """
    unique: dict[tuple, SomalierFlag] = {}
    for sf in sg_flags:
        for flag in sf.flags:
            unique[flag_identity(flag, flag_sg_key(flag, sf.sg_id))] = flag

    all_flags = list(unique.values())
    active = [f for f in all_flags if not f.resolved]
    active_by_category = {key: sum(1 for f in active if CATEGORY_KEYS.get(f.category) == key) for key in CATEGORY_ORDER}

    return {
        'total_sgs': total_sgs,
        'active_flags': len(active),
        'active_by_category': active_by_category,
        'families_affected': families_affected,
        'resolved_flags': sum(1 for f in all_flags if f.resolved),
    }


def category_chips(groups: list[FamilyGroup]) -> list[dict]:
    """Filter chips for the categories actually present, counting families rather than flags."""
    chips = []
    for key in CATEGORY_ORDER:
        count = sum(1 for group in groups if group.counts.get(key))
        if count:
            chips.append({'key': key, 'label': CATEGORY_LABELS[key], 'count': count})
    return chips


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def render_report(
    dataset: str,
    active_groups: list[FamilyGroup],
    resolved_groups: list[FamilyGroup],
    summary: dict,
    generated_at: str | None = None,
) -> str:
    """Render the Somalier flags report HTML using the Jinja template."""
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
        autoescape=jinja2.select_autoescape(['html', 'xml']),
    )
    template = env.get_template('somalier_flags_overview.html.jinja')
    chips = category_chips(active_groups)
    return template.render(
        dataset=dataset,
        generated_at=generated_at or datetime.now(tz=UTC).isoformat(timespec='seconds'),
        active_groups=active_groups,
        resolved_groups=resolved_groups,
        summary=summary,
        category_chips=chips,
        inline_flag_limit=INLINE_FLAG_LIMIT,
        show_filter_bar=len(chips) > 1 or len(active_groups) > MIN_GROUPS_FOR_FILTER_BAR,
        cross_family_mark=CROSS_FAMILY_MARK,
        dash=DASH,
    )


def construct_summary_message(
    dataset: str,
    *,
    flags_html_url: str,
    somalier_html_url: str,
    seq_type: str,
    seq_tech: str,
    summary: dict,
    previous_analysis: dict | None,
):
    """Construct a Slack message with a concise summary and a link to the report."""


def main(dataset: str, output_html: str, base_output_html: str, flags_html_url: str, somalier_html_url: str) -> None:
    """Query Metamist for Somalier flags and generate a Somalier flags HTML report."""

    dataset = dataset_for_access_level(dataset)
    seq_type = config_retrieve(['workflow', 'sequencing_type'])
    seq_tech = config_retrieve(['workflow', 'sequencing_technology'])

    logging_prefix = f'{dataset} ({seq_type} | {seq_tech})'

    logger.info(f'{logging_prefix} :: Querying Metamist for Somalier flags')
    started = perf_counter()
    response = query(DATASET_SGS_QUERY, variables={'dataset': dataset, 'seqType': seq_type, 'seqTech': seq_tech})
    sequencing_groups = response['project']['sequencingGroups']
    logger.info(
        f'{logging_prefix} :: Found {len(sequencing_groups)} sequencing groups in {perf_counter() - started:.1f}s'
    )

    sg_flags = collect_somalier_flags(sequencing_groups)
    flagged = [sf for sf in sg_flags if sf.flags]
    logger.info(f'{logging_prefix} :: {len(flagged)} sequencing groups have Somalier flags.')

    # Includes the far member of every pairwise flag, which is not otherwise in the flagged list.
    infos = get_sg_infos(referenced_sg_ids(flagged))

    groups = group_by_family(flagged, infos)
    active_groups, resolved_groups = split_active_resolved(groups, infos)
    summary = summarise_flags(flagged, total_sgs=len(sequencing_groups), families_affected=len(active_groups))

    logger.info(
        f'{logging_prefix} :: Rendering {summary["active_flags"]} active flag(s) across '
        f'{len(active_groups)} famil{"y" if len(active_groups) == 1 else "ies"}'
    )
    started = perf_counter()
    html = render_report(dataset, active_groups, resolved_groups, summary=summary)
    logger.info(f'{logging_prefix} :: Rendered report in {perf_counter() - started:.1f}s')

    with to_path(base_output_html).open('w') as f:
        f.write(html)
    logger.info(f'{logging_prefix} :: Wrote Somalier flags report to {base_output_html}')

    with to_path(output_html).open('w') as f:
        f.write(html)
    logger.info(f'{logging_prefix} :: Wrote timestamped Somalier flags report to {output_html}')

    # Register results in Metamist manually to capture all dataset SGs in scope, not just the input_cohorts SGs
    meta = {
        'stage': STAGE_NAME,
        'dataset': dataset,
        'sequencing_type': seq_type,
        'sequencing_technology': seq_tech,
        'summary': summary,
    }

    create_new(
        project=dataset,
        output=output_html,
        analysis_type='web',
        sgs=[sg['id'] for sg in sequencing_groups],
        meta=meta,
    )
    logger.info(f'{logging_prefix} :: Registered web analysis for {len(sequencing_groups)} SG(s)')
    logger.info(f'{logging_prefix} :: {flags_html_url}')
    logger.info(f'{logging_prefix} :: {somalier_html_url}')

    meta.pop('summary')
    construct_summary_message(
        dataset,
        flags_html_url=flags_html_url,
        somalier_html_url=somalier_html_url,
        seq_type=seq_type,
        seq_tech=seq_tech,
        summary=summary,
        previous_analysis=get_previous_analysis(dataset, meta),
    )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Metamist dataset/project name')
    parser.add_argument('--output-html', required=True, help='gs:// path to HTML (namespaced by timestamp)')
    parser.add_argument('--base-output-html', required=True, help='gs:// path to HTML (fixed, not namespaced)')
    parser.add_argument('--flags-html-url', required=True, help='Clickable URL for the Somalier flags report')
    parser.add_argument('--somalier-html-url', required=True, help='Clickable URL for the original Somalier report')
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        output_html=args.output_html,
        base_output_html=args.base_output_html,
        flags_html_url=args.flags_html_url,
        somalier_html_url=args.somalier_html_url,
    )
