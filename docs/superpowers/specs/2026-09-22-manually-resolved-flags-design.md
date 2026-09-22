# Manually resolved flags: design

Status: approved 2026-09-22. Builds on `2026-09-14-somalier-relatedness-report-design.md` for the flag data model and `2026-09-18-relatedness-report-final-touches-design.md` for the report's sections.

## Purpose

A curator who reviews a flag and decides it is a known, accepted finding currently has no way to say so. Every run rewrites `SG.meta['somalier_flags']` from scratch, and reconciliation is mechanical: a flag whose identity key reappears is written back with `resolved=False` (`record_somalier_flags.py:272-277`). Editing the record by hand in Metamist buys one run of quiet before the pipeline undoes it.

The cases this is for are the ones where nothing is going to change. A pedigree that is known to be wrong and will not be corrected, a relatedness conflict already investigated and explained, a sex inference mismatch on a sample nobody is going to re-sequence. These keep leading the report, and their presence at the top makes it harder to see the findings that are new.

What this adds is a small CLI that records a manual resolution against a specific flag, a reconciler that honours it across runs, and a report that leaves those flags out.

## What a manual resolution is

Three new fields on the `SomalierFlag` base dataclass (`utils.py:134`), inherited by all three categories:

```python
manually_resolved: bool = False
manual_resolution_reason: str | None = None
manual_resolution_by: str | None = None
```

No separate date field. A manual resolution sets the existing `resolved=True` and `resolution_date` (`utils.py:145-146`) exactly as an automatic one does, and `manually_resolved` is the only thing distinguishing the two. Every field defaults, so flags written before this change deserialise unchanged and pick the fields up the next time they are rewritten, the same way `sequencing_group_key` was introduced.

The marker goes **on the flag record itself**, not into a separate table of suppressions keyed by SG pair. This is what makes the binding strict for free: reconciliation already matches on the full identity key, so a pair whose measured relationship changes produces a different record, and that record has no marker on it.

## Identifying a flag

`(sequencing_group_key, category)` is enough, and no flag ID or hash is needed.

`sequencing_group_key` is the sorted, underscore-joined SG IDs the flag involves (`record_somalier_flags.py:43-55`), so `CPG001_CPG002` for a pairwise flag and `CPG001` for a per-SG one. Both pairwise categories are recorded against the sorted-first SG of the pair (`check_pedigree.py:269-281`, `check_self_relatedness.py:190-194`), which means the key also names the SG whose meta holds the flag: it is the first element. Passing an SG ID separately would be redundant.

The pair is unique among **unresolved** flags, which is all the CLI ever needs to act on. Two properties give that:

- One flag per key per category per run. `check_pedigree.py:233` iterates one row per pair out of `pairs_df`, and `check_self_relatedness.py:189` one per failing pair, so a run cannot emit two flags sharing a key and category.
- When a pair's measurement changes, the old identity key goes absent from the new run and is marked resolved (`record_somalier_flags.py:250-256`) in the same pass that adds the new one. So the older record is never left unresolved alongside its replacement.

It is *not* unique across the whole stored list, because resolution history accumulates. A pair flagged in March as `parent-child` and now as `unrelated` are two entries with the same key and category, one resolved and one not. The CLI narrows to unresolved flags for that reason, and refuses to act if it still somehow finds more than one rather than picking one.

## The CLI

`src/rd_qc/scripts/resolve_somalier_flag.py`, with a `resolve_somalier_flag` entry point in `pyproject.toml` alongside `run_workflow`. It runs locally, against the curator's own Metamist credentials, not inside a Hail Batch job.

```bash
resolve_somalier_flag --dataset my-dataset \
  --sg-ids CPG001 CPG002 \
  --category relatedness_mismatch \
  --reason "pedigree known wrong, family declined correction" \
  --reviewer ef
```

`--sg-ids` takes one or two IDs in any order and sorts them into the key, so the curator can read the two IDs straight off a report row without worrying which one owns the flag. `--category` is one of the three category strings, passed as `choices` so a typo fails at parse time rather than matching nothing.

It reads the owning SG's meta, narrows to unresolved flags matching the key and category, prints the matched flag in full, and writes the whole `somalier_flags` list back on `y`. `--yes` skips the prompt for scripted use.

It exits non-zero without writing when no unresolved flag matches, when more than one does, or when the match is already manually resolved. In the first two cases it prints the flags it did find for that key, including resolved ones, since "I resolved that last month" and "I have the wrong SG" look identical from an empty result otherwise.

`--unresolve` is the inverse. It clears the three new fields along with `resolved` and `resolution_date`, so the flag is active again and returns to the report on the next run. Reopening a mistake should not need a hand-edit in Metamist. `--reason` and `--reviewer` are required on both paths and rejected if empty or whitespace, because an unexplained suppression is worse than no suppression.

## Shared access to the meta key

`DATASET_SG_META_QUERY` and `SG_META_MUTATION` live in `record_somalier_flags.py:13-40` today. Both move to a new `src/rd_qc/flag_store.py`, which exposes a read and a write over `meta['somalier_flags']` for the reconciler and the CLI to share.

The mutation overwrites the whole list, so a second copy written independently is a second chance to drop every flag on an SG. Worth keeping in one place, and `flag_store.py` also gives the CLI's tests one thing to patch.

## Reconciliation

The behaviour change is one branch in each of the three reconcile loops (`record_somalier_flags.py:122`, `:188`, `:249`), ahead of the `compare_*` call:

```python
if flag_key not in new_flags:        # unchanged: absent, so resolve or leave resolved
elif flag.get('manually_resolved'):  # new: refresh measured values, keep the resolution
elif compare_...(flag, new_flag):    # unchanged: same unresolved issue, refresh values
else:                                # unchanged: differs or reappeared, overwrite
```

Without it a manually-resolved flag whose key reappears falls through `compare_*`, which requires `not resolved` (`record_somalier_flags.py:66`, `:81`, `:98`), into the `else` branch, where `flag.update(new_flag)` sets `resolved=False` and leaves the three manual fields sitting on an active flag. That record would be self-contradictory and the suppression would last exactly one run.

The measured-value refresh is currently written out three times with a different key list each time (`:137-142`, `:204-206`, `:265-269`), and the new branch needs the same thing. Rather than adding a fourth copy, that becomes one `refresh_measured_values(flag, new_flag, keys)` helper called from both branches. The three reconcile functions are otherwise left alone: they are near-identical and could plausibly collapse into one generic function, but that is a bigger change than this feature needs and it would make the diff hard to review.

`manually_resolved` joins the stats dict so the per-SG log line reports how many were held.

### When the finding genuinely goes away

A manually-resolved flag whose key stops appearing (someone corrected the pedigree after all) keeps its marker and stays hidden, by the existing first branch: `if not flag['resolved']` is already false, so it is left as-is.

The alternative considered was clearing the marker at that point, so the flag would drop into the report's resolved-backlog section as an ordinary past incident. Rejected for being more machinery in exchange for a section nobody is waiting on. The consequence accepted with it: an accepted finding that later resolves on its own never shows up in the backlog, and the only record of it is in Metamist.

## Report

One skip at the top of the loop in `collect_somalier_flags` (`somalier_flags_report.py:447`), logged at debug level.

That function is the only place flags enter the report, so a flag dropped there is absent from all four sections, from every count in the summary cards, and from the Slack post, with no other file touched. Manually resolved flags are invisible in the report, not shown as a count.

A count in the summary bar was the other option. It would mean threading a number through `summarise_flags`, `render_report`, the template, `_headline_lines` and the previous-run delta comparison (`somalier_flags_report.py:1039`), which is a lot of plumbing for one number. It stays available as a later addition if the invisibility turns out to be a problem, since the records themselves are untouched in Metamist.

## Testing

New `test/test_resolve_somalier_flag.py`, patching `flag_store` the way `test_record_somalier_flags.py` already patches the query:

- a matching unresolved flag is found and the written payload carries the three fields plus `resolved` and `resolution_date`
- no match, two unresolved matches, and an already-manually-resolved match each exit non-zero and write nothing
- `--unresolve` clears all five fields
- an empty or whitespace `--reason` or `--reviewer` is rejected

Added to `test_record_somalier_flags.py`, covering the lifecycle across runs:

- the same identity key reappearing leaves `resolved` and the manual fields intact, with measured values refreshed
- a changed identity key for the same pair is added as a new active flag, and the manually-resolved record stays resolved
- a legacy flag dict with none of the new keys still reconciles

Added to `test_somalier_flags_report.py`: a manually-resolved flag appears in neither the active nor the resolved sections, and does not move any summary count.

## Out of scope

Bulk or batch resolution. One flag per invocation, which keeps the confirmation prompt meaningful.

Any change to which flags get produced. `check_pedigree.py` and `check_self_relatedness.py` are untouched, so a suppressed finding is still measured and still recorded, and the Slack messages those checks send are unaffected.

Surfacing reasons or reviewers anywhere. The report hides these flags rather than annotating them, so the reason text stays in Metamist and out of a shared HTML file.
