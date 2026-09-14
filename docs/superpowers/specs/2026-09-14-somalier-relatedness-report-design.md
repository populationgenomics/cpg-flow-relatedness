# Somalier relatedness report: design

Status: approved 2026-09-14. Supersedes the structural recommendation in `cpg-flow-align-genotype/docs/somalier-relatedness-report-handoff.md` section 7b.

## Purpose

`src/rd_qc/scripts/somalier_flags_report.py` queries every sequencing group in a dataset, reads the Somalier flags recorded in `SG.meta['somalier_flags']`, and renders a single self-contained HTML dashboard. The flag-producing and flag-recording steps already exist (`check_pedigree.py`, `check_self_relatedness.py`, `record_somalier_flags.py`). This spec covers only the report, which is the third and final step.

The report is for collaborators looking at relatedness problems in their own data, so it leads with family and participant identifiers and demotes CPG sequencing group IDs to small muted monospace.

## Current state

Both the script and the template were forked from the SG QC report in `cpg-flow-align-genotype` and left partly adapted. Neither runs today.

`src/rd_qc/templates/somalier_flags_overview.html.jinja` is a verbatim copy of `sg_qc_overview.html.jinja`. It reads `f.source`, `f.severity`, `f.severity_label`, `f.metric_label`, `f.section_label`, `sg.count_summary`, `summary.active_fail`, `summary.active_cram`, `generated_at`, `unresolved`, `resolved` and `active_metrics`. `render_report` (`somalier_flags_report.py:254`) passes only `dataset`, `reports` and `summary`, so every one of those resolves to undefined and the page renders as an empty shell.

Four defects need fixing before anything else works:

1. `somalier_flags_report.py:228` calls `SomalierFlag(**flag)` on every flag. `SomalierFlag` (`utils.py:113`) accepts only `category`, `date`, `ar_guid`, `resolved` and `resolution_date`, so any real flag carrying `provided=` or `sg_id_1=` raises `TypeError`.
2. `somalier_flags_report.py:240` computes `active_relatedness_flags` by iterating `all_flags` rather than `active`, so resolved pedigree flags are counted as active. The two sibling counts on lines 238 and 239 iterate `active` correctly, which is what makes this look like a slip rather than intent.
3. `render_report` never passes `generated_at`, which the template reads at line 242.
4. `src/rd_qc/jobs/somalier_flags_report.py:28` invokes `python3 -m align_genotype.scripts.somalier_flags_report`. The module lives at `rd_qc.scripts.somalier_flags_report`, so the Hail Batch job fails immediately.
5. `config.dataset_for_access_level` was added in cpg-utils 5.7.2, and the only constraint on it was the loose transitive one from `cpg-flow~=1.3`. This was a local-only breakage, not a CI one: CI resolves fresh via `pip install .[test]` and was green at 0509628, but the untracked local `uv.lock` had pinned 5.7.0, which broke the report's import along with all four `stages.py` call sites (lines 25, 47, 75, 102) and three `test_stage_scoping.py` tests. Pinned `cpg-utils>=5.7.2` in `pyproject.toml` regardless, because the floor is real and neither CI nor the Docker build reads a lockfile, so pyproject is the only place that can enforce it.

## Two findings that shape the design

### Pairwise flags are not duplicated across sequencing groups

The handoff document's main warning (section 7b) was that a pairwise flag gets written into both members' meta, so naive per-SG rendering double-counts it and a dedup pass is required. That does not hold in this repo. Both producers deliberately record a pair flag against the lexicographically first sequencing group only:

- `check_self_relatedness.py:132-140` sorts the pair into `s1, s2` then appends to `flags_by_sg_id[s1]`, with the comment "Only necessary to register the flag for the first SG in each pair".
- `check_pedigree.py:186-199` does the same, sorting `s1, s2` and appending to `relatedness_flags_by_sg_id[s1]`.

So each pairwise flag exists exactly once across the dataset's meta and no dedup pass is needed on read. What is needed instead is metadata for the partner: `sg_id_2` is never fetched today, so the second sample's family and participant IDs are unavailable to the renderer.

### No severity is recorded on any flag

None of the three dataclasses in `utils.py` carries a severity field, so the copied template's whole fail/warn colour system has nothing to bind to. `check_pedigree.py:213-219` does distinguish the two directions, treating expected-related-inferred-unrelated as an error and expected-unrelated-inferred-related (above a `relatedness > 0.1` cutoff) as a warning, but that distinction only ever reaches the Slack message text. It is not persisted.

Decision: drop severity from the report entirely rather than deriving it at render time or adding it at the source. Every active flag renders identically with a single red left border. Deriving a fail/warn split in the renderer would invent a clinical judgement the pipeline never made, and adding the field to the dataclasses would need a fallback for every flag already sitting in Metamist without it. Neither is worth it for a colour.

## Page model

Group by family. A family is the natural review unit for relatedness, and it is the only grouping under which a pedigree mismatch and a sex mismatch on the same individual appear next to each other.

Structure is two sections, mirroring the QC report:

- `⚠ Families with flags`, the headline, groups rendered open so flags are visible without clicking.
- `✓ Resolved, past incidents`, the same macro at 72% opacity.

In both sections the flag lines themselves are always visible, one line each, and it is the detail row (measured values plus per-SG metadata) that starts collapsed. Hiding the flag lines in the resolved section would leave a backlog you cannot read without clicking every row, which defeats the point of showing it.

A family with both active and resolved flags appears in both sections, showing only the relevant flags in each. Colour semantics carry over from the QC report: red only when active flags exist, green only for the all-clear banner, grey for resolved and neutral counts. A flag count is never green.

Header cards are `SGs scanned`, `Active flags`, `Families affected`, `Resolved`.

## Data model

Three frozen dataclasses in `somalier_flags_report.py`, all display-layer only.

```python
@dataclass(frozen=True)
class FlagRow:
    category: str            # 'sex_inference_mismatch' | 'self_relatedness_mismatch' | 'relatedness_mismatch'
    category_key: str        # 'sex' | 'self' | 'pedigree', for filter chips and data-* attributes
    category_label: str      # 'Sex inference' | 'Self-relatedness' | 'Pedigree relatedness'
    identity: tuple          # dataset-wide dedup and count key
    sg_key: str              # the resolved sequencing_group_key, used to collect a group's SGInfos
    subject: str             # 'PID_C / CPG004' or 'CPG004 <-> CPG005'
    subject_detail: str      # participant behind the pair: 'PID_B' or 'PID_C <-> PID_D'
    result: str              # 'provided female / inferred male'
    details: tuple[tuple[str, str], ...]   # metric label to formatted value, for the expanded table
    cross_family: str | None # the other family, set only on duplicated cross-family rows
    resolved: bool
    date_short: str
    date_full: str
    resolution_date_short: str
    resolution_date_full: str
    search_blob: str         # lowercased, rolled up into the group's own search text


@dataclass(frozen=True)
class FamilyGroup:
    key: str                       # 'FAM02' | 'participant:PID_B' | 'sg:CPG009'
    label: str                     # 'FAM02' | '(no family) / PID_B'
    flags: tuple[FlagRow, ...]
    sg_infos: tuple[SGInfo, ...]   # every SG referenced by this group's flags
    counts: dict[str, int]         # keyed by category_key
    count_summary: str             # '3 flags (1 sex, 2 pedigree)'
    search_blob: str               # lowercased family, participant, sample and SG ids
```

`SGInfo` (`somalier_flags_report.py:86`) stays as it is. `SGReport` is removed, since the family group replaces it as the render unit.

`identity` is what keeps the counts honest. A cross-family row appears in two `FamilyGroup`s but carries one `identity`, so `summarise_flags` deduplicates on it and reports the flag once:

| category | identity |
| --- | --- |
| sex | `('sex', sequencing_group_key, provided, inferred)` |
| self | `('self', sequencing_group_key, participant_external_id, threshold)` |
| pedigree | `('pedigree', sequencing_group_key, expected_relationship, inferred_relationship)` |

These mirror the recording script's own identity keys (`record_somalier_flags.py:97`, `:161`, `:222`), so the report and the resolved/unresolved lifecycle agree on what counts as the same flag. The measured values (`relatedness`, `ibs0`, `ibs2`, and the sex statistics) are excluded from identity in both places because they drift between relate runs.

## The persisted sequencing group key

`SomalierFlag` gains one field, which is the flag's answer to "which sequencing groups am I about?":

```python
sequencing_group_key: str = ''
```

The value is the sorted, underscore-joined set of SG IDs the flag involves, so `CPG001` for a sex flag and `CPG002_CPG003` for either kind of pairwise flag. `sg_ids_tag` (`utils.py:348`) already computes exactly this and is reused rather than reimplemented.

It is populated in `record_somalier_flags.py`, not in the two producer scripts. That file already rebuilds the per-category identity tuples three times over (lines 97, 161, 222), so it is where identity logic belongs, and setting the key during reconciliation means one code path covers flags from both producers.

The field defaults to `''` so that every flag already sitting in Metamist still deserialises. The report therefore needs a fallback for legacy flags: when `sequencing_group_key` is empty, derive it from `sg_id_1` and `sg_id_2` for pairwise categories, or from the owning SG ID for sex flags. That fallback is the only place the report looks at those fields directly.

This closes both of the TODOs added on 2026-09-14, at `utils.py:117` and `somalier_flags_report.py:300`. The second one is closed because `referenced_sg_ids` can now split the key on `_` to learn every SG a flag touches, rather than needing per-category knowledge of where the partner ID lives.

## Pipeline

```text
main
 |- query DATASET_SGS_QUERY                     unchanged
 |- collect_somalier_flags  -> list[SgFlags]    rewrite: dispatch on category
 |- referenced_sg_ids       -> set[str]         new: split sequencing_group_key on '_'
 |- get_sg_infos            -> dict[str, SGInfo] unchanged, already keyed by id
 |- build_flag_rows         -> list[FlagRow]    new: one branch per category
 |- group_by_family         -> list[FamilyGroup] new: the core of this design
 |- split_active_resolved   -> (active, resolved)
 |- summarise_flags         -> dict             rewrite: deduplicate on identity
 |- render_report
```

`collect_somalier_flags` dispatches through a module-level map instead of the current bare `SomalierFlag(**flag)`:

```python
FLAG_CLASSES = {
    'sex_inference_mismatch': SomalierSexInferenceFlag,
    'self_relatedness_mismatch': SomalierSelfRelatednessFlag,
    'relatedness_mismatch': SomalierRelatednessFlag,
}
```

An unrecognised `category`, or a flag whose fields do not fit its class, logs a warning and is skipped. A report that dies on one malformed meta entry is worse than one that renders the other ninety-nine.

`referenced_sg_ids` must union the owning SG with `sg_id_1` and `sg_id_2` from every pairwise flag, because the partner is not otherwise in the query list and its family is needed for grouping.

### Group assignment

One rule governs the whole page.

| category | primary key | fallback chain |
| --- | --- | --- |
| sex | family of the owning SG, from Metamist | `participant:<pid>`, then `sg:<sg_id>` |
| self | family of `sg_id_1`, from Metamist | `participant:<participant_external_id>` |
| pedigree | families of `sg_id_1` and `sg_id_2`, from Metamist, falling back to the flag's `family_external_id` | `participant:`, then `sg:`, per member |

For pedigree flags, when the two resolved families differ and both are known, the row is emitted into both groups with `cross_family` set to the other family's label and rendered with a `↗ cross-family` marker. The duplication is deliberate. A pair that is expected unrelated but inferred related across two families is the classic cross-family sample swap, and both families' reviewers need to see it. The `identity` key stops it being counted twice.

Metamist is preferred over the flag's stored `family_external_id` because that field is already lossy. `check_pedigree.py:194` collapses it to `fam1 or fam2 or 'unknown'`, which discards the fact that the pair straddled two families at all. Metamist is the only place the second family survives.

`split_active_resolved` takes the single list of groups and returns two lists of freshly built `FamilyGroup`s: one holding only each group's unresolved flags, one holding only its resolved flags, with `counts`, `count_summary` and `sg_infos` recomputed against the filtered flag set. A group with nothing left after filtering is dropped from that side. This is what lets a family appear in both sections showing only the flags relevant to each.

Within each section, groups sort by flag count descending, then by label, so the worst families surface first. Flags within a group sort by category, then subject.

## Template

Same file, `src/rd_qc/templates/somalier_flags_overview.html.jinja`. The CSS shell, header cards, filter bar JavaScript and `id-*` identifier styling are reused as they are. Changes:

- Delete the severity system: `.badge-fail`, `.badge-warn`, `.sev-*`, `.sev-count-*`, the severity chips, the severity branch in `applyFilters()`, and the Failing and Warning header cards.
- Delete the source system: `.badge-source`, the source chips, `summary.active_cram` and `summary.active_gvcf`. There is no CRAM/GVCF axis in relatedness flags.
- The outer table row becomes a family group: label, `count_summary`, and every flag line inline (subject, result, and the cross-family marker where set). Expanding a group reveals the full per-flag metric table plus per-SG metadata for each SG involved.
- Filter chips become the three `category_key` values with per-family counts. Search covers `search_blob`, so it matches either member of a pair.
- Render the filter bar only when there is something to filter, meaning more than one category or more than five groups.

The `{% macro %}` must be defined before it is called, so it stays immediately after `<body>`. Autoescape is on, so display strings are built in Python rather than with `|safe`.

Number formatting keeps `_fmt_num` from the QC report: integers stay integers, values at or above 1 get two decimal places, values below 1 get two significant figures, trailing zeros are stripped. So `0.616722` renders as `0.62` and `4821.0` as `4821`. Dates are truncated to `YYYY-MM-DD` for display with the full ISO timestamp in a hover `title`.

Category labels come from a small hand-maintained map with a fallback to the raw key, following the same reasoning as `METRIC_LABELS` in the QC report. The raw category strings are not collaborator-friendly and friendly names are not available anywhere upstream.

## Mock data and offline render

Verify the look offline rather than against the live API, which dodges the transient Metamist GraphQL stalls described in the handoff document's section 6.

`test/fixtures/somalier_flags.py` holds hand-built `meta['somalier_flags']` blobs in exact database shape, covering:

- an active sex mismatch, and a resolved one
- a sex flag with `provided='unknown'`
- an active self-relatedness flag
- a pedigree mismatch within one family
- a pedigree mismatch across two families, which should appear in both groups
- a family carrying both active and resolved flags
- a sequencing group with no family at all, exercising the participant fallback

`testing_scripts/render_mock_somalier_report.py` feeds those through `collect_somalier_flags`, `group_by_family` and `render_report`, then writes `somalier_flags_report.html` plus a second all-clear render into `--output-dir` (the system temp directory by default, rather than a hardcoded `/tmp`, which trips ruff's `S108`). It also prints the resolved group structure, so the page can be sanity checked without opening a browser. No Metamist, fully deterministic.

## Tests

`test/test_somalier_flags_report.py`, pure functions only, no database. The `test/**/*.py` per-file-ignores block in `pyproject.toml` exempted only `S101`, so `PLR2004` was added: an expected count is the whole point of an assertion. `RUF001` turned out not to be needed, since ruff does not treat `↔` or `·` as ambiguous.

Cases:

- category dispatch produces the right dataclass for each of the three categories
- an unknown category is skipped and logged, not raised
- a cross-family pedigree flag lands in both family groups
- that same flag is counted once in `summarise_flags`
- a self-relatedness flag groups under the family resolved from Metamist
- a self-relatedness flag with no family falls back to a participant group
- active and resolved split correctly, with a family appearing in both sections
- `_fmt_num` across integer, above-one and below-one inputs
- rendered HTML contains the expected identifiers and result strings
- the all-clear banner renders when there are no active flags

Commands:

```bash
uv run --with pytest python -m pytest test/ -q
uv run --with ruff ruff check src/ test/ testing_scripts/
```

## Out of scope

`construct_summary_message` (`somalier_flags_report.py:264`) is an empty stub called from `main`. It returns `None` harmlessly so the report works, but the Slack summary does nothing and `get_previous_analysis` is queried for no reason. Left alone here.

Adding a *severity* field to the flag dataclasses stays out of scope, per the decision above. Note that this is narrower than it was: `sequencing_group_key` is now in scope, so `utils.py` and `record_somalier_flags.py` do both get touched, just not for severity.
