# Manually resolved flags implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a curator mark a Somalier flag as a reviewed, accepted finding so it stays resolved across pipeline runs and drops out of the HTML report.

**Architecture:** A manual resolution is three fields written onto the flag record itself in `SG.meta['somalier_flags']`, not a separate suppression list. A local CLI writes them against one flag, reconciliation grows one branch that holds such a flag resolved when its finding recurs, and the report skips them where it reads flags out of meta. Because the marker lives on the record and reconciliation matches on the full identity key, a changed measurement is a different record and surfaces as a new active flag with no marker on it.

**Tech Stack:** Python 3.10+, dataclasses, `metamist.graphql` (`gql`/`query`), loguru, argparse, pytest, ruff (line length 120, single quotes).

**Spec:** `docs/superpowers/specs/2026-09-22-manually-resolved-flags-design.md`. Read it first.

---

## Background an engineer needs

Somalier QC findings are stored as a list of dicts under `meta['somalier_flags']` on each sequencing group (SG) in Metamist. There are three categories, each with its own dataclass in `src/rd_qc/utils.py`: `sex_inference_mismatch`, `self_relatedness_mismatch`, `relatedness_mismatch`.

Every pipeline run rewrites that whole list per SG (`record_somalier_flags.py`). Reconciliation compares stored flags against the run's new flags on an **identity key** per category, then marks absent ones resolved, refreshes recurring ones, and appends new ones. Identity keys deliberately exclude the measured values (`relatedness`, `ibs0`, `ibs2`, depth ratios) because those drift between runs.

Pairwise flags (a pair of SGs) are recorded against the **sorted-first** SG of the pair only, and carry a `sequencing_group_key` field holding the sorted, underscore-joined IDs (`CPG001_CPG002`). Per-SG flags key on the one SG. So the key also tells you which SG's meta holds the flag: the first element.

Nothing here touches Hail Batch. The CLI runs on a laptop; the reconciler runs inside a Batch job via `python3 -m rd_qc.scripts.record_somalier_flags` (`jobs/relate.py:214`).

**Run tests with:** `uv run pytest test -q` from the repo root. `pythonpath = ['src']` is set in `pyproject.toml`, so no install step is needed. Ruff runs via pre-commit, and is not a project dependency, so invoke it with the version pre-commit pins: `uv run --with ruff==0.15.19 ruff check src test` and `uv run --with ruff==0.15.19 ruff format src test`.

## File structure

| File | Responsibility |
|---|---|
| `src/rd_qc/utils.py` (modify, ~line 134) | Three new fields on the `SomalierFlag` base dataclass, inherited by all three categories. |
| `src/rd_qc/flag_store.py` (create) | The only module that reads or writes `meta['somalier_flags']`: the two GraphQL operations, `sequencing_group_key`, and read/write helpers. |
| `src/rd_qc/scripts/record_somalier_flags.py` (modify) | Reconciliation. Loses its GraphQL operations to `flag_store`, gains the manual-resolution branch and a shared measured-value refresh helper. |
| `src/rd_qc/scripts/resolve_somalier_flag.py` (create) | The curator-facing CLI. Pure selection and transformation functions, plus a thin `main` that reads, confirms, and writes. |
| `src/rd_qc/scripts/somalier_flags_report.py` (modify, ~line 447) | One skip in `collect_somalier_flags`, the single place flags enter the report. |
| `test/test_flag_store.py` (create) | Read and write helpers against a patched `query`. |
| `test/test_resolve_somalier_flag.py` (create) | CLI selection, refusal cases, transformations, and `main` exit codes. |
| `test/test_record_somalier_flags.py` (modify) | Fixture repoint plus the manual-resolution lifecycle across runs. |
| `test/test_somalier_flags_report.py` (modify) | Held flags reach no section and no count. |
| `pyproject.toml`, `README.md` (modify) | Console script entry point, docs, version bump. |

Tasks 1 through 4 are the feature. Task 2 is a pure refactor with no behaviour change and must leave the suite green on its own.

---

### Task 1: Manual-resolution fields on the flag dataclass

**Files:**
- Modify: `src/rd_qc/utils.py:134-146`
- Test: `test/test_record_somalier_flags.py`

- [ ] **Step 1: Write the failing tests**

Add to the imports at the top of `test/test_record_somalier_flags.py`:

```python
from rd_qc.utils import SomalierRelatednessFlag
```

Append to the end of the file:

```python
def test_manual_resolution_fields_default_to_absent():
    """A flag nobody has reviewed carries the fields, unset, so every record has the same shape."""
    flag = SomalierRelatednessFlag(
        category='relatedness_mismatch',
        sg_id_1='CPG1',
        sg_id_2='CPG2',
        family_external_id='FAM1',
        expected_relationship='siblings',
        inferred_relationship='unrelated',
        relatedness=0.02,
        ibs0=900,
        ibs2=100,
    )

    assert flag.manually_resolved is False
    assert flag.manual_resolution_reason is None
    assert flag.manual_resolution_by is None


def test_a_flag_stored_before_the_manual_fields_existed_still_deserialises():
    """The fixture dict has none of the new keys, which is what Metamist holds for older flags."""
    flag = SomalierRelatednessFlag(**relatedness_flag())

    assert flag.manually_resolved is False
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_record_somalier_flags.py -q -k manual`
Expected: FAIL, `AttributeError: 'SomalierRelatednessFlag' object has no attribute 'manually_resolved'`

- [ ] **Step 3: Add the fields**

In `src/rd_qc/utils.py`, inside `class SomalierFlag`, immediately after the existing `resolution_date` field (line 146):

```python
    resolved: bool = False
    resolution_date: str | None = None
    # A resolution a curator recorded by hand with resolve_somalier_flag, rather than one
    # reconciliation inferred from the finding going away. Held across runs even while the finding
    # recurs, and skipped by the report. `resolved` and `resolution_date` are set alongside these,
    # so a reader that only knows about automatic resolution still sees a resolved flag.
    manually_resolved: bool = False
    manual_resolution_reason: str | None = None
    manual_resolution_by: str | None = None
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test -q`
Expected: PASS, whole suite green. The fields are keyword-only with defaults, so nothing that builds a flag today needs changing.

- [ ] **Step 5: Commit**

```bash
git add src/rd_qc/utils.py test/test_record_somalier_flags.py
git commit -m "feat: add manual resolution fields to the somalier flag base class"
```

---

### Task 2: Extract `flag_store.py`

Pure refactor. No behaviour changes, no new fields used. The point is that the whole-list mutation exists once, because a second independent copy of it is a second chance to wipe every flag on an SG.

**Files:**
- Create: `src/rd_qc/flag_store.py`
- Modify: `src/rd_qc/scripts/record_somalier_flags.py:1-55`, `:376-384`, `:421-434`
- Test: `test/test_flag_store.py` (create), `test/test_record_somalier_flags.py:56-70`

- [ ] **Step 1: Write the failing tests**

Create `test/test_flag_store.py`:

```python
"""
Tests for the one module that reads and writes SG meta's `somalier_flags` list.

The mutation overwrites the whole list, so both the reconciler and the resolve CLI go through here.
Metamist is patched out: what matters is the shape of what would be sent.
"""

from rd_qc import flag_store

FLAG = {'category': 'relatedness_mismatch', 'sg_id_1': 'CPG1', 'sg_id_2': 'CPG2', 'resolved': False}


def fake_dataset(monkeypatch, sequencing_groups: list[dict]) -> dict:
    """Patch `query` to answer the dataset read, and capture any mutation variables."""
    captured: dict = {}

    def fake_query(_query, variables=None) -> dict:
        captured.update(variables or {})
        return {'project': {'sequencingGroups': sequencing_groups}}

    monkeypatch.setattr(flag_store, 'query', fake_query)
    return captured


def test_read_sg_flags_returns_the_stored_list(monkeypatch):
    fake_dataset(monkeypatch, [{'id': 'CPG1', 'meta': {'somalier_flags': [FLAG]}}])

    assert flag_store.read_sg_flags('my-dataset', 'CPG1') == [FLAG]


def test_read_sg_flags_returns_empty_for_a_sequencing_group_with_no_flags(monkeypatch):
    """
    An SG that exists but has never been flagged is not the same as a missing SG.

    Its meta is populated with other things, so this is the common real shape: the key is absent
    rather than the meta being empty.
    """
    fake_dataset(monkeypatch, [{'id': 'CPG1', 'meta': {'sequencing_type': 'genome'}}])

    assert flag_store.read_sg_flags('my-dataset', 'CPG1') == []


def test_read_sg_flags_returns_empty_for_a_sequencing_group_with_null_meta(monkeypatch):
    """Metamist returns meta as null rather than {} for some SGs, which must not raise."""
    fake_dataset(monkeypatch, [{'id': 'CPG1', 'meta': None}])

    assert flag_store.read_sg_flags('my-dataset', 'CPG1') == []


def test_read_sg_flags_returns_none_for_a_sequencing_group_not_in_the_dataset(monkeypatch):
    """Distinguishable from the empty case so a caller can say 'wrong dataset' rather than 'no flags'."""
    fake_dataset(monkeypatch, [{'id': 'CPG1', 'meta': {}}])

    assert flag_store.read_sg_flags('my-dataset', 'CPG9') is None


def test_write_sg_flags_sends_the_whole_list_under_the_meta_key(monkeypatch):
    captured = fake_dataset(monkeypatch, [])

    flag_store.write_sg_flags('my-dataset', 'CPG1', [FLAG])

    assert captured == {'dataset': 'my-dataset', 'sgId': 'CPG1', 'sgMeta': {'somalier_flags': [FLAG]}}


def test_sequencing_group_key_joins_a_pair_in_sorted_order():
    assert flag_store.sequencing_group_key({'sg_id_1': 'CPG2', 'sg_id_2': 'CPG1'}, 'CPG2') == 'CPG1_CPG2'


def test_sequencing_group_key_falls_back_to_the_owning_sg_for_a_per_sg_flag():
    assert flag_store.sequencing_group_key({'provided': 'M', 'inferred': 'F'}, 'CPG1') == 'CPG1'
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_flag_store.py -q`
Expected: FAIL at collection, `ModuleNotFoundError: No module named 'rd_qc.flag_store'`

- [ ] **Step 3: Create the module**

Create `src/rd_qc/flag_store.py`:

```python
"""
Reading and writing the `somalier_flags` list on a sequencing group's meta in Metamist.

Both the pipeline's reconciler and the resolve CLI write this list, and the mutation replaces it
wholesale rather than patching entries, so a second independently-written copy of it is a second
chance to drop every flag on an SG. It lives here once.
"""

from rd_qc.utils import sg_ids_tag

from metamist.graphql import gql, query

SOMALIER_FLAGS_KEY = 'somalier_flags'

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


def sequencing_group_key(flag: dict, sg_id: str) -> str:
    """
    Sorted, underscore-joined SG IDs that this flag involves.

    Pairwise flags (self-relatedness, relatedness) are about two SGs but are recorded against
    only the first of the pair, so this key is what lets a reader work out which SGs a flag
    touches without needing per-category knowledge of where the partner ID lives. Per-SG flags
    (sex inference) key on the SG that owns them.
    """
    sg_id_1, sg_id_2 = flag.get('sg_id_1'), flag.get('sg_id_2')
    if sg_id_1 and sg_id_2:
        return sg_ids_tag([sg_id_1, sg_id_2])
    return sg_id


def read_dataset_sg_meta(dataset: str) -> list[dict]:
    """Every sequencing group in the dataset, as `{'id': ..., 'meta': ...}`."""
    response = query(DATASET_SG_META_QUERY, variables={'dataset': dataset})
    return response['project']['sequencingGroups']


def read_sg_flags(dataset: str, sg_id: str) -> list[dict] | None:
    """
    The stored Somalier flags for one sequencing group.

    `None` means the dataset has no such SG, which is a different problem from an SG that has
    never been flagged, and the two deserve different error messages.
    """
    for sg in read_dataset_sg_meta(dataset):
        if sg['id'] == sg_id:
            return (sg.get('meta') or {}).get(SOMALIER_FLAGS_KEY, [])
    return None


def write_sg_flags(dataset: str, sg_id: str, flags: list[dict]) -> None:
    """Replace one sequencing group's whole Somalier flag list."""
    query(
        SG_META_MUTATION,
        variables={'dataset': dataset, 'sgId': sg_id, 'sgMeta': {SOMALIER_FLAGS_KEY: flags}},
    )
```

- [ ] **Step 4: Point the reconciler at it**

In `src/rd_qc/scripts/record_somalier_flags.py`, delete `DATASET_SG_META_QUERY` and `SG_META_MUTATION` (lines 13-40) and the `sequencing_group_key` function (lines 43-55), delete the `from metamist.graphql import gql, query` import, and replace the `rd_qc.utils` import block with:

```python
from rd_qc.flag_store import read_dataset_sg_meta, sequencing_group_key, write_sg_flags
from rd_qc.utils import SomalierRelatednessFlag, SomalierSelfRelatednessFlag, SomalierSexInferenceFlag
```

Replace the mutation call at the end of `reconcile_sg_somalier_flags` (lines 376-384):

```python
    write_sg_flags(dataset, sg_id, [asdict(flag) for flag in final_flags])
```

And the dataset read in `main` (lines 421-423):

```python
    # Query the sequencing groups for the given dataset
    sequencing_groups = read_dataset_sg_meta(dataset)
```

Delete the now-unused local variable `somalier_flags_key = 'somalier_flags'` at line 304 and use the imported constant in the one place it is still read:

```python
    current_somalier_flags: list[dict] = (sg['meta'] or {}).get(SOMALIER_FLAGS_KEY, [])
```

adding `SOMALIER_FLAGS_KEY` to the `rd_qc.flag_store` import line.

- [ ] **Step 5: Repoint the existing test fixture**

The `written_meta` fixture patches `record_somalier_flags.query`, which no longer exists there. In `test/test_record_somalier_flags.py`, add to the imports:

```python
from rd_qc import flag_store
```

and change the one line inside the fixture (line 69) from `monkeypatch.setattr(record_somalier_flags, 'query', fake_query)` to:

```python
    monkeypatch.setattr(flag_store, 'query', fake_query)
```

The captured variables keep the same shape, so every existing assertion on `written_meta['sgMeta']['somalier_flags']` is unaffected.

- [ ] **Step 6: Run the whole suite**

Run: `uv run pytest test -q`
Expected: PASS. This task changes no behaviour, so a failure here is a refactor mistake, not a feature gap.

- [ ] **Step 7: Commit**

```bash
git add src/rd_qc/flag_store.py src/rd_qc/scripts/record_somalier_flags.py test/test_flag_store.py test/test_record_somalier_flags.py
git commit -m "refactor: move somalier flag meta access into flag_store"
```

---

### Task 3: Reconciliation holds a manual resolution

**Files:**
- Modify: `src/rd_qc/scripts/record_somalier_flags.py` (all three reconcile loops)
- Test: `test/test_record_somalier_flags.py`

- [ ] **Step 1: Write the failing tests**

Append to `test/test_record_somalier_flags.py`:

```python
def manually_resolved(flag: dict, **overrides: object) -> dict:
    """`flag` as the resolve CLI leaves it: resolved, with the reviewer's reason attached."""
    return (
        flag
        | {
            'resolved': True,
            'resolution_date': RESOLVED_EARLIER,
            'manually_resolved': True,
            'manual_resolution_reason': 'pedigree known wrong',
            'manual_resolution_by': 'ef',
        }
        | overrides
    )


def test_manually_resolved_flag_stays_resolved_when_the_finding_recurs(written_meta):
    """
    The whole point of the feature: a run that measures the same thing again must not reopen it.

    Without the manual branch this flag falls through compare_* into the overwrite branch, which
    sets resolved=False and leaves the manual fields on an active flag.
    """
    held = manually_resolved(relatedness_flag())
    recurrence = relatedness_flag(date=TODAY, relatedness=0.05, ibs0=850, ibs2=120)

    reconcile(current_flags=[held], new_flags=[recurrence])

    written = flags_by_category(written_meta)['relatedness_mismatch']
    assert written['resolved'] is True
    assert written['manually_resolved'] is True
    assert written['manual_resolution_by'] == 'ef'
    assert written['manual_resolution_reason'] == 'pedigree known wrong'
    assert written['resolution_date'] == RESOLVED_EARLIER, 'the reviewer resolved it, not this run'
    assert written['date'] == FIRST_SEEN, 'a held issue keeps its first-detected date'
    assert (written['relatedness'], written['ibs0'], written['ibs2']) == (0.05, 850, 120)


def test_a_changed_finding_is_not_suppressed_by_a_manual_resolution(written_meta):
    """
    Binding is strict: the marker is on one record, so a different measurement surfaces unheld.

    This is the safety property. A pair accepted as parent-child must not stay quiet when the
    genotypes start saying unrelated.
    """
    held = manually_resolved(relatedness_flag())
    reinferred = relatedness_flag(inferred_relationship='parent-child')

    reconcile(current_flags=[held], new_flags=[reinferred])

    written = written_meta['sgMeta']['somalier_flags']
    assert len(written) == 2

    by_inferred = {flag['inferred_relationship']: flag for flag in written}
    assert by_inferred['unrelated']['manually_resolved'] is True
    assert by_inferred['parent-child']['resolved'] is False
    assert by_inferred['parent-child']['manually_resolved'] is False


def test_a_manually_resolved_flag_stays_held_once_the_finding_disappears(written_meta):
    """An accepted finding that later goes away keeps its marker and stays out of the report."""
    held = manually_resolved(relatedness_flag())

    reconcile(current_flags=[held], new_flags=[sex_flag()])

    written = flags_by_category(written_meta)['relatedness_mismatch']
    assert written['manually_resolved'] is True
    assert written['resolution_date'] == RESOLVED_EARLIER, 'the resolution date is not re-stamped'


def test_manual_resolution_is_held_for_every_category(written_meta):
    """The branch is duplicated across three reconcilers, so all three get pinned."""
    held_sex = manually_resolved(sex_flag())

    reconcile(current_flags=[held_sex], new_flags=[sex_flag(mean_depth=29.0)])

    written = flags_by_category(written_meta)['sex_inference_mismatch']
    assert written['resolved'] is True
    assert written['manually_resolved'] is True
    assert written['mean_depth'] == 29.0, 'measured values still refresh while held'
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_record_somalier_flags.py -q -k "held or suppressed or recurs or every_category"`
Expected: FAIL. `test_manually_resolved_flag_stays_resolved_when_the_finding_recurs` fails on `assert written['resolved'] is True` because the overwrite branch reopened it.

- [ ] **Step 3: Add the refresh helper and the constants**

In `src/rd_qc/scripts/record_somalier_flags.py`, directly below the imports:

```python
# The measured fields each category refreshes on a flag it is keeping. Identity fields are absent
# by construction: a change in any of those makes a different flag, not an update to this one.
SEX_MEASURED_FIELDS = ('mean_depth', 'x_het_ratio', 'x_depth_ratio', 'y_depth_ratio', 'p_middling_ab')
SELF_RELATEDNESS_MEASURED_FIELDS = ('relatedness', 'ibs0', 'ibs2')
# `verdict` is derived from the measured values, so it refreshes with them.
RELATEDNESS_MEASURED_FIELDS = ('relatedness', 'ibs0', 'ibs2', 'verdict')


def refresh_measured_values(flag: dict, new_flag: dict, fields: tuple[str, ...]) -> None:
    """
    Copy this run's measurements onto a flag being kept, leaving identity and resolution alone.

    Mutates in place, matching how the reconcilers below already update these dicts. A field the
    new flag does not carry is left as it was rather than blanked.
    """
    flag.update({field: new_flag[field] for field in fields if field in new_flag})
```

- [ ] **Step 4: Add the manual branch to the sex reconciler**

In `reconcile_sg_somalier_sex_inference_flags`, change the stats line to:

```python
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0, 'held': 0}
```

Then insert a new branch between the `if flag_key not in ...` block and the `elif compare_...` block, and rewrite the retained branch's body to use the helper, so the chain reads:

```python
        elif flag.get('manually_resolved'):
            # Reviewed and accepted by a curator. The finding is still here, so take this run's
            # measurements but leave the resolution alone. Without this branch the flag falls
            # through compare_* (which requires an unresolved flag) into the overwrite below,
            # which would reopen it every run.
            refresh_measured_values(flag, new_somalier_sex_inference_flags_by_key[flag_key], SEX_MEASURED_FIELDS)
            logger.info(
                f"{sg_id} :: {report} flag '{flag['provided']}-{flag['inferred']}' "
                f'manually resolved by {flag.get("manual_resolution_by")}; held.'
            )
            stats['held'] += 1
        elif compare_somalier_sex_inference_flag(flag, new_somalier_sex_inference_flags_by_key[flag_key]):
            # Same unresolved issue is still present: refresh the measured value and but keep resolution status.
            # Identity (provided/inferred) is unchanged so this counts as 'retained', not 'updated'.
            refresh_measured_values(flag, new_somalier_sex_inference_flags_by_key[flag_key], SEX_MEASURED_FIELDS)
            logger.info(
                f"{sg_id} :: {report} flag '{flag['provided']}-{flag['inferred']}' "
                'remains unresolved (value refreshed).'
            )
            stats['retained'] += 1
```

- [ ] **Step 5: Add the same branch to the self-relatedness reconciler**

In `reconcile_sg_somalier_self_relatedness_flags`, the stats line gains `'held': 0` as above, and the chain becomes:

```python
        elif flag.get('manually_resolved'):
            # See the sex reconciler above: held across runs, measurements still refreshed.
            refresh_measured_values(
                flag,
                new_somalier_self_relatedness_flags_by_key[flag_key],
                SELF_RELATEDNESS_MEASURED_FIELDS,
            )
            logger.info(
                f"{sg_id} :: {report} flag '{flag['sg_id_1']}-{flag['sg_id_2']}' "
                f'manually resolved by {flag.get("manual_resolution_by")}; held.'
            )
            stats['held'] += 1
        elif compare_somalier_self_relatedness_flag(flag, new_somalier_self_relatedness_flags_by_key[flag_key]):
            # Same unresolved issue is still present: refresh the measured value and
            # but keep resolution status. Identity (sg_id_1/sg_id_2/participant_external_id/threshold)
            # is unchanged so this counts as 'retained', not 'updated'.
            refresh_measured_values(
                flag,
                new_somalier_self_relatedness_flags_by_key[flag_key],
                SELF_RELATEDNESS_MEASURED_FIELDS,
            )
            logger.info(
                f"{sg_id} :: {report} flag '{flag['sg_id_1']}-{flag['sg_id_2']}' remains unresolved (value refreshed)."
            )
            stats['retained'] += 1
```

- [ ] **Step 6: Add the same branch to the relatedness reconciler**

In `reconcile_sg_somalier_relatedness_flags`, the stats line gains `'held': 0`, and the chain becomes:

```python
        elif flag.get('manually_resolved'):
            # See the sex reconciler above: held across runs, measurements still refreshed.
            refresh_measured_values(flag, new_somalier_relatedness_flags_by_key[flag_key], RELATEDNESS_MEASURED_FIELDS)
            logger.info(
                f"{sg_id} :: {report} flag '{flag['category']}' "
                f'manually resolved by {flag.get("manual_resolution_by")}; held.'
            )
            stats['held'] += 1
        elif compare_somalier_relatedness_flag(flag, new_somalier_relatedness_flags_by_key[flag_key]):
            # Same unresolved issue is still present: refresh the measured value and but keep resolution status.
            # Identity (sg_id_1/sg_id_2/family_external_id/expected_relationship/inferred_relationship)  # noqa: ERA001
            # is unchanged so this counts as 'retained', not 'updated'.
            refresh_measured_values(flag, new_somalier_relatedness_flags_by_key[flag_key], RELATEDNESS_MEASURED_FIELDS)
            logger.info(f"{sg_id} :: {report} flag '{flag['category']}' remains unresolved (value refreshed).")
            stats['retained'] += 1
```

- [ ] **Step 7: Carry `held` through the aggregate and the log line**

In `reconcile_sg_somalier_flags`, the aggregate stats initialiser (line 342) becomes:

```python
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0, 'held': 0}
```

The three `stats = {k: stats[k] + ....get(k, 0) for k in stats}` lines need no change, since they iterate the aggregate's keys. Extend the closing log line so a held flag is visible in the job log:

```python
    logger.info(
        f'{sg_id} :: Recorded {len(final_flags)} {report} flags in Metamist. '
        f'Resolved: {stats["resolved"]}, Retained: {stats["retained"]}, '
        f'Updated: {stats["updated"]}, Added: {stats["added"]}, '
        f'Manually resolved: {stats["held"]}'
    )
```

- [ ] **Step 8: Run the tests to verify they pass**

Run: `uv run pytest test -q`
Expected: PASS, including the four pre-existing reconciliation tests. `test_recurring_flag_is_retained_with_refreshed_measurements` is the one that proves the helper did not change the retained path.

- [ ] **Step 9: Commit**

```bash
git add src/rd_qc/scripts/record_somalier_flags.py test/test_record_somalier_flags.py
git commit -m "feat: hold manually resolved somalier flags across runs"
```

---

### Task 4: The report skips held flags

**Files:**
- Modify: `src/rd_qc/scripts/somalier_flags_report.py:443-457`
- Test: `test/test_somalier_flags_report.py`

- [ ] **Step 1: Write the failing tests**

Add to `test/test_somalier_flags_report.py`, after the existing `test_collect_handles_sg_with_no_meta_at_all`:

```python
# A manual resolution as the resolve CLI writes it, for spreading over a fixture flag.
HELD = {
    'resolved': True,
    'resolution_date': RESOLVED_ON,
    'manually_resolved': True,
    'manual_resolution_reason': 'pedigree known wrong',
    'manual_resolution_by': 'ef',
}


def test_collect_skips_a_manually_resolved_flag():
    groups = [{'id': 'CPG001', 'meta': {'somalier_flags': [sex_flag('CPG001', 'M', 'F') | HELD]}}]

    collected = collect_somalier_flags(groups)

    assert collected[0].flags == ()


def test_manually_resolved_flags_reach_no_section_and_no_count():
    """
    Held flags are invisible, not merely de-emphasised.

    Not a conflict, not a refinement, not resolved history, and not in any summary number, so a
    reader cannot mistake an accepted finding for a fixed one.
    """
    _, _, _, baseline = run_pipeline()
    assert baseline['active_flags'] > 0, 'the baseline fixture must have something to hide'

    held_groups = [
        {
            'id': sg['id'],
            'meta': {'somalier_flags': [flag | HELD for flag in (sg['meta'] or {}).get('somalier_flags', [])]},
        }
        for sg in MOCK_SEQUENCING_GROUPS
    ]

    flagged, active, resolved, summary = run_pipeline(held_groups)

    assert flagged == []
    assert (active, resolved) == ([], [])
    assert summary['active_flags'] == 0
    assert summary['active_conflicts'] == 0
    assert summary['active_refinements'] == 0
    assert summary['resolved_flags'] == 0
```

`RESOLVED_ON` already exists in `test/fixtures/somalier_flags.py`; add it to the `from fixtures.somalier_flags import (...)` block at the top of the test file.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_somalier_flags_report.py -q -k "manually_resolved"`
Expected: FAIL. `test_collect_skips_a_manually_resolved_flag` fails because `collected[0].flags` has one `SomalierSexInferenceFlag` in it.

- [ ] **Step 3: Add the skip**

In `collect_somalier_flags`, immediately inside the `for raw in meta.get('somalier_flags') or []:` loop, before the category dispatch:

```python
        for raw in meta.get('somalier_flags') or []:
            if (raw or {}).get('manually_resolved'):
                # Reviewed and accepted by a curator, so it is not a finding this report is for.
                # Dropped here rather than downstream because this is the only way flags enter the
                # report: skipping here keeps it out of all four sections and every summary count.
                logger.debug(f'{sg["id"]} :: skipping manually resolved Somalier flag')
                continue
            category = (raw or {}).get('category') or ''
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test -q`
Expected: PASS. Watch `test_summary_counts` and `test_all_clear_banner_when_nothing_is_active` in particular: the fixtures carry no held flags, so their numbers must not move.

- [ ] **Step 5: Commit**

```bash
git add src/rd_qc/scripts/somalier_flags_report.py test/test_somalier_flags_report.py
git commit -m "feat: keep manually resolved flags out of the somalier report"
```

---

### Task 5: CLI selection and transformation

The pure half of the CLI, tested without any Metamist or terminal involvement. `main` arrives in Task 6.

**Files:**
- Create: `src/rd_qc/scripts/resolve_somalier_flag.py`
- Test: `test/test_resolve_somalier_flag.py`

- [ ] **Step 1: Write the failing tests**

Create `test/test_resolve_somalier_flag.py`:

```python
"""
Tests for the manual flag resolution CLI.

Selection is the part worth testing. Resolution history means one key and category can match
several stored flags, so the CLI has to pick the single unresolved one or refuse, and it must never
guess. Metamist and the terminal are patched out.
"""

from rd_qc.scripts import resolve_somalier_flag as cli

FIRST_SEEN = '2026-01-01T00:00:00+00:00'
RESOLVED_EARLIER = '2026-02-02T00:00:00+00:00'
NOW = '2026-09-22T00:00:00+00:00'
REASON = 'pedigree known wrong, family declined correction'
REVIEWER = 'ef'


def relatedness_flag(**overrides: object) -> dict:
    """A relatedness_mismatch flag as stored in SG meta, owned by CPG1 and keyed on the pair."""
    return {
        'category': 'relatedness_mismatch',
        'sequencing_group_key': 'CPG1_CPG2',
        'date': FIRST_SEEN,
        'resolved': False,
        'resolution_date': None,
        'manually_resolved': False,
        'manual_resolution_reason': None,
        'manual_resolution_by': None,
        'sg_id_1': 'CPG1',
        'sg_id_2': 'CPG2',
        'family_external_id': 'FAM1',
        'expected_relationship': 'siblings',
        'inferred_relationship': 'unrelated',
        'relatedness': 0.02,
        'ibs0': 900,
        'ibs2': 100,
    } | overrides


def held(**overrides: object) -> dict:
    """A flag already manually resolved."""
    return relatedness_flag(
        resolved=True,
        resolution_date=RESOLVED_EARLIER,
        manually_resolved=True,
        manual_resolution_reason=REASON,
        manual_resolution_by=REVIEWER,
        **overrides,
    )


def select(flags: list[dict], *, unresolve: bool = False) -> tuple[dict | None, str]:
    return cli.select_target(flags, 'CPG1_CPG2', 'relatedness_mismatch', 'CPG1', unresolve=unresolve)


def test_flag_key_is_order_independent():
    """A curator reads two IDs off a report row and should not have to know which sorts first."""
    assert cli.flag_key_of(['CPG2', 'CPG1']) == 'CPG1_CPG2'
    assert cli.flag_key_of(['CPG1', 'CPG2']) == 'CPG1_CPG2'
    assert cli.flag_key_of(['CPG1']) == 'CPG1'


def test_owning_sg_is_the_first_of_the_key():
    """Pairwise flags are recorded against the sorted-first SG, so the key names its own owner."""
    assert cli.owning_sg_id('CPG1_CPG2') == 'CPG1'
    assert cli.owning_sg_id('CPG1') == 'CPG1'


def test_the_one_unresolved_flag_is_selected():
    target = relatedness_flag()

    selected, problem = select([target])

    assert selected is target
    assert problem == ''


def test_resolution_history_does_not_make_the_match_ambiguous():
    """The same pair flagged differently in the past is a separate, already-resolved record."""
    history = relatedness_flag(inferred_relationship='parent-child', resolved=True, resolution_date=RESOLVED_EARLIER)
    target = relatedness_flag()

    selected, _ = select([history, target])

    assert selected is target


def test_a_flag_for_another_pair_is_not_selected():
    other = relatedness_flag(sequencing_group_key='CPG1_CPG3', sg_id_2='CPG3')

    selected, problem = select([other])

    assert selected is None
    assert 'No relatedness_mismatch flag stored' in problem


def test_a_flag_of_another_category_is_not_selected():
    sex = {'category': 'sex_inference_mismatch', 'sequencing_group_key': 'CPG1', 'resolved': False}

    selected, problem = select([sex])

    assert selected is None
    assert 'No relatedness_mismatch flag stored' in problem


def test_a_legacy_flag_without_a_recorded_key_still_matches():
    """Flags written before `sequencing_group_key` existed are matched on their SG IDs instead."""
    legacy = relatedness_flag()
    del legacy['sequencing_group_key']

    selected, _ = select([legacy])

    assert selected is legacy


def test_an_already_held_flag_is_refused_with_the_reviewer_and_reason():
    selected, problem = select([held()])

    assert selected is None
    assert 'already manually resolved' in problem
    assert REVIEWER in problem
    assert '--unresolve' in problem


def test_an_auto_resolved_flag_is_refused_and_the_stored_flags_are_shown():
    """'I resolved that last month' and 'I have the wrong SG' must not read the same."""
    selected, problem = select([relatedness_flag(resolved=True, resolution_date=RESOLVED_EARLIER)])

    assert selected is None
    assert 'No unresolved relatedness_mismatch flag' in problem
    assert 'parent-child' not in problem
    assert 'siblings' in problem, 'the stored flag is printed back so the curator can see it'


def test_two_unresolved_candidates_are_refused_rather_than_guessed_at():
    """
    Reconciliation should make this unreachable, since a changed key resolves the old flag in the
    same pass. If it ever happens, picking one silently would hide a live finding.
    """
    selected, problem = select([relatedness_flag(), relatedness_flag(expected_relationship='parent-child')])

    assert selected is None
    assert 'refusing to guess' in problem


def test_unresolve_selects_the_held_flag_not_the_active_one():
    target = held()
    active = relatedness_flag(inferred_relationship='parent-child')

    selected, _ = select([target, active], unresolve=True)

    assert selected is target


def test_unresolve_is_refused_when_nothing_is_held():
    selected, problem = select([relatedness_flag()], unresolve=True)

    assert selected is None
    assert 'No manually resolved relatedness_mismatch flag' in problem


def test_with_manual_resolution_sets_the_five_fields_and_touches_nothing_else():
    original = relatedness_flag()

    updated = cli.with_manual_resolution(original, REASON, REVIEWER, NOW)

    assert updated['resolved'] is True
    assert updated['resolution_date'] == NOW
    assert updated['manually_resolved'] is True
    assert updated['manual_resolution_reason'] == REASON
    assert updated['manual_resolution_by'] == REVIEWER
    assert updated['date'] == FIRST_SEEN, 'the first-detected date is not the resolution date'
    assert updated['inferred_relationship'] == 'unrelated', 'identity is untouched'
    assert original['resolved'] is False, 'the stored flag is not mutated'


def test_without_manual_resolution_clears_every_trace():
    updated = cli.without_manual_resolution(held())

    assert updated['resolved'] is False
    assert updated['resolution_date'] is None
    assert updated['manually_resolved'] is False
    assert updated['manual_resolution_reason'] is None
    assert updated['manual_resolution_by'] is None


def test_replace_flag_swaps_one_entry_and_keeps_the_order():
    first = relatedness_flag(inferred_relationship='parent-child', resolved=True)
    target = relatedness_flag()
    replacement = cli.with_manual_resolution(target, REASON, REVIEWER, NOW)

    written = cli.replace_flag([first, target], target, replacement)

    assert written == [first, replacement]


def test_replace_flag_matches_on_identity_not_equality():
    """Two history records can be equal dicts; only the selected one may be rewritten."""
    twin = relatedness_flag()
    target = relatedness_flag()
    replacement = cli.with_manual_resolution(target, REASON, REVIEWER, NOW)

    written = cli.replace_flag([twin, target], target, replacement)

    assert written[0] is twin
    assert written[1] is replacement
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_resolve_somalier_flag.py -q`
Expected: FAIL at collection, `ModuleNotFoundError: No module named 'rd_qc.scripts.resolve_somalier_flag'`

- [ ] **Step 3: Write the module**

Create `src/rd_qc/scripts/resolve_somalier_flag.py`:

```python
"""
Record or clear a manual resolution on a single Somalier flag in Metamist.

A finding a curator has reviewed and accepted (a pedigree known to be wrong, a mismatch already
investigated and explained) reappears in the report every run, because reconciliation reopens any
flag it measures again. Marking it manually resolved holds it resolved across runs and keeps it out
of the report, without stopping the checks from measuring it.

Run locally, against your own Metamist credentials:

    resolve_somalier_flag --dataset my-dataset --sg-ids CPG001 CPG002 \
        --category relatedness_mismatch --reason "pedigree known wrong" --reviewer ef

Take the SG IDs off the report row in either order. Add --unresolve to reopen a flag.
"""

import json
from argparse import ArgumentParser
from datetime import UTC, datetime

from loguru import logger

from rd_qc.flag_store import read_sg_flags, sequencing_group_key, write_sg_flags
from rd_qc.utils import sg_ids_tag

CATEGORIES = ('sex_inference_mismatch', 'self_relatedness_mismatch', 'relatedness_mismatch')

EXIT_OK = 0
EXIT_REFUSED = 1
EXIT_BAD_ARGS = 2


def flag_key_of(sg_ids: list[str]) -> str:
    """The `sequencing_group_key` for these SG IDs, given in either order."""
    return sg_ids_tag(sg_ids)


def owning_sg_id(sg_key: str) -> str:
    """
    The sequencing group whose meta holds a flag with this key.

    Pairwise flags are recorded against the sorted-first SG of the pair and the key is the sorted
    join, so the owner is the first element.
    """
    return sg_key.split('_')[0]


def describe(flags: list[dict]) -> str:
    """One indented JSON line per flag, for printing candidates back to the curator."""
    return '\n'.join(f'  {json.dumps(flag, sort_keys=True)}' for flag in flags)


def matching_flags(flags: list[dict], sg_key: str, category: str, owner: str) -> list[dict]:
    """
    Every stored flag for this key and category, resolved or not.

    Falls back to deriving the key for flags recorded before `sequencing_group_key` existed, the
    same way the report does.
    """
    return [
        flag
        for flag in flags
        if flag.get('category') == category
        and (flag.get('sequencing_group_key') or sequencing_group_key(flag, owner)) == sg_key
    ]


def select_target(
    flags: list[dict],
    sg_key: str,
    category: str,
    owner: str,
    *,
    unresolve: bool,
) -> tuple[dict | None, str]:
    """
    The single flag to act on, or `None` and the reason why not.

    Resolving acts on an unresolved flag, unresolving on a manually resolved one. Anything other
    than exactly one candidate is refused rather than guessed at: resolution history means a key
    and category can match several stored records, and picking one silently could hide a live
    finding. Refusals print the stored flags, so 'I resolved that last month' and 'I have the wrong
    sequencing group' do not read the same.
    """
    stored = matching_flags(flags, sg_key, category, owner)
    if not stored:
        return None, f'No {category} flag stored for {sg_key}. Check the SG IDs and the category.'

    if unresolve:
        candidates = [flag for flag in stored if flag.get('manually_resolved')]
        if not candidates:
            return None, f'No manually resolved {category} flag for {sg_key}. Stored flags:\n{describe(stored)}'
    else:
        candidates = [flag for flag in stored if not flag.get('resolved', False)]
        if not candidates:
            already_held = [flag for flag in stored if flag.get('manually_resolved')]
            if already_held:
                flag = already_held[0]
                return None, (
                    f'{category} for {sg_key} is already manually resolved, on '
                    f'{flag.get("resolution_date")} by {flag.get("manual_resolution_by")}: '
                    f'{flag.get("manual_resolution_reason")}. Use --unresolve to reopen it.'
                )
            return None, f'No unresolved {category} flag for {sg_key}. Stored flags:\n{describe(stored)}'

    if len(candidates) > 1:
        return None, (
            f'{len(candidates)} candidate {category} flags for {sg_key}; refusing to guess.\n{describe(candidates)}'
        )
    return candidates[0], ''


def with_manual_resolution(flag: dict, reason: str, reviewer: str, now: str) -> dict:
    """
    A copy of `flag`, resolved by hand.

    `resolved` and `resolution_date` are set alongside the manual fields so a reader that knows
    nothing about manual resolution still sees a resolved flag. Identity fields and the
    first-detected `date` are untouched.
    """
    return flag | {
        'resolved': True,
        'resolution_date': now,
        'manually_resolved': True,
        'manual_resolution_reason': reason,
        'manual_resolution_by': reviewer,
    }


def without_manual_resolution(flag: dict) -> dict:
    """A copy of `flag`, active again, with every trace of the manual resolution cleared."""
    return flag | {
        'resolved': False,
        'resolution_date': None,
        'manually_resolved': False,
        'manual_resolution_reason': None,
        'manual_resolution_by': None,
    }


def replace_flag(flags: list[dict], target: dict, replacement: dict) -> list[dict]:
    """
    A new list with `target` swapped for `replacement`, matched on identity rather than equality.

    Two history records for one pair can be equal dicts, and only the selected one may be rewritten.
    """
    return [replacement if flag is target else flag for flag in flags]
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test/test_resolve_somalier_flag.py -q`
Expected: PASS, all 16 tests.

- [ ] **Step 5: Commit**

```bash
git add src/rd_qc/scripts/resolve_somalier_flag.py test/test_resolve_somalier_flag.py
git commit -m "feat: add flag selection logic for manual somalier resolutions"
```

---

### Task 6: CLI entry point

**Files:**
- Modify: `src/rd_qc/scripts/resolve_somalier_flag.py`
- Test: `test/test_resolve_somalier_flag.py`

- [ ] **Step 1: Write the failing tests**

Add `import pytest` to the top of `test/test_resolve_somalier_flag.py`, above the `from rd_qc.scripts import ...` line (ruff's E402 rejects a mid-file import). Then append the rest:

```python
@pytest.fixture
def metamist(monkeypatch):
    """
    Patch the flag store, confirm the prompt, and expose what would have been written.

    `written` stays empty when nothing was written, which is what every refusal must produce.
    """
    state: dict = {'stored': [], 'written': {}}

    def fake_read(_dataset: str, _sg_id: str) -> list[dict]:
        return state['stored']

    def fake_write(dataset: str, sg_id: str, flags: list[dict]) -> None:
        state['written'] = {'dataset': dataset, 'sg_id': sg_id, 'flags': flags}

    monkeypatch.setattr(cli, 'read_sg_flags', fake_read)
    monkeypatch.setattr(cli, 'write_sg_flags', fake_write)
    monkeypatch.setattr(cli, 'confirmed', lambda: True)
    return state


def run(metamist_state: dict, **overrides: object) -> int:
    """Invoke main with the usual arguments, overriding as needed."""
    kwargs = {
        'dataset': 'my-dataset',
        'sg_ids': ['CPG2', 'CPG1'],
        'category': 'relatedness_mismatch',
        'reason': REASON,
        'reviewer': REVIEWER,
        'unresolve': False,
        'assume_yes': False,
    } | overrides
    return cli.main(**kwargs)


def test_resolving_writes_the_marked_flag_against_the_owning_sg(metamist):
    metamist['stored'] = [relatedness_flag()]

    assert run(metamist) == 0

    written = metamist['written']
    assert written['dataset'] == 'my-dataset'
    assert written['sg_id'] == 'CPG1', 'pairwise flags live on the sorted-first SG'
    assert len(written['flags']) == 1
    assert written['flags'][0]['manually_resolved'] is True
    assert written['flags'][0]['manual_resolution_by'] == REVIEWER
    assert written['flags'][0]['manual_resolution_reason'] == REASON
    assert written['flags'][0]['resolution_date'], 'a resolution date is stamped'


def test_the_rest_of_the_list_is_written_back_untouched(metamist):
    """The mutation replaces the whole list, so anything dropped here is lost from Metamist."""
    history = relatedness_flag(inferred_relationship='parent-child', resolved=True)
    metamist['stored'] = [history, relatedness_flag()]

    assert run(metamist) == 0

    assert metamist['written']['flags'][0] == history


def test_unresolving_clears_the_marker(metamist):
    metamist['stored'] = [held()]

    assert run(metamist, unresolve=True) == 0

    written = metamist['written']['flags'][0]
    assert written['manually_resolved'] is False
    assert written['resolved'] is False


def test_a_refused_selection_writes_nothing(metamist):
    metamist['stored'] = [held()]

    assert run(metamist) == 1
    assert metamist['written'] == {}


def test_an_sg_missing_from_the_dataset_writes_nothing(metamist, monkeypatch):
    monkeypatch.setattr(cli, 'read_sg_flags', lambda *_: None)

    assert run(metamist) == 1
    assert metamist['written'] == {}


def test_declining_the_prompt_writes_nothing(metamist, monkeypatch):
    metamist['stored'] = [relatedness_flag()]
    monkeypatch.setattr(cli, 'confirmed', lambda: False)

    assert run(metamist) == 1
    assert metamist['written'] == {}


def test_assume_yes_skips_the_prompt(metamist, monkeypatch):
    metamist['stored'] = [relatedness_flag()]

    def refuse() -> bool:
        raise AssertionError('the prompt must not be reached with --yes')

    monkeypatch.setattr(cli, 'confirmed', refuse)

    assert run(metamist, assume_yes=True) == 0


@pytest.mark.parametrize('blank', ['', '   '])
def test_an_empty_reason_or_reviewer_is_rejected_before_anything_is_read(metamist, blank):
    """An unexplained suppression is worse than none, so this fails at the door."""
    metamist['stored'] = [relatedness_flag()]

    assert run(metamist, reason=blank) == 2
    assert run(metamist, reviewer=blank) == 2
    assert metamist['written'] == {}


def test_the_reason_and_reviewer_are_stored_stripped(metamist):
    metamist['stored'] = [relatedness_flag()]

    run(metamist, reason=f'  {REASON}  ', reviewer='  ef  ')

    written = metamist['written']['flags'][0]
    assert written['manual_resolution_reason'] == REASON
    assert written['manual_resolution_by'] == 'ef'
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest test/test_resolve_somalier_flag.py -q -k "writes or unresolving or prompt or rejected or stripped"`
Expected: FAIL, `AttributeError: module 'rd_qc.scripts.resolve_somalier_flag' has no attribute 'main'`

- [ ] **Step 3: Add `confirmed`, `main` and the argument parser**

Append to `src/rd_qc/scripts/resolve_somalier_flag.py`:

```python
def confirmed() -> bool:
    """Ask before writing. Anything other than an explicit 'y' is a no."""
    return input('Write this change? [y/N] ').strip().lower() == 'y'


def main(
    dataset: str,
    sg_ids: list[str],
    category: str,
    reason: str,
    reviewer: str,
    *,
    unresolve: bool = False,
    assume_yes: bool = False,
) -> int:
    """
    Resolve or reopen one flag. Returns the process exit code and writes at most once.

    Nothing is written unless exactly one flag matched and the change was confirmed.
    """
    if not reason.strip() or not reviewer.strip():
        logger.error('--reason and --reviewer must both be non-empty.')
        return EXIT_BAD_ARGS

    sg_key = flag_key_of(sg_ids)
    owner = owning_sg_id(sg_ids)

    flags = read_sg_flags(dataset, owner)
    if flags is None:
        logger.error(f'{owner} is not a sequencing group in {dataset}.')
        return EXIT_REFUSED

    target, problem = select_target(flags, sg_key, category, owner, unresolve=unresolve)
    if target is None:
        logger.error(problem)
        return EXIT_REFUSED

    action = 'Reopening' if unresolve else 'Manually resolving'
    logger.info(f'{action} this flag on {owner} in {dataset}:\n{describe([target], sg_key)}')
    if not assume_yes and not confirmed():
        logger.info('Aborted, nothing written.')
        return EXIT_REFUSED

    now = datetime.now(tz=UTC).isoformat(timespec='seconds')
    replacement = (
        without_manual_resolution(target)
        if unresolve
        else with_manual_resolution(target, reason.strip(), reviewer.strip(), now)
    )
    write_sg_flags(dataset, owner, replace_flag(flags, target, replacement))
    logger.info(f'Wrote {len(flags)} flags back to {owner} in {dataset}.')
    return EXIT_OK


def cli_main() -> int:
    parser = ArgumentParser(description='Record or clear a manual resolution on one Somalier flag.')
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument(
        '--sg-ids',
        nargs='+',
        required=True,
        help='The SG IDs the flag involves, in either order: one for a sex flag, two for a pair',
    )
    parser.add_argument('--category', required=True, choices=CATEGORIES, help='Flag category')
    parser.add_argument('--reason', required=True, help='Why this finding is accepted as-is')
    parser.add_argument('--reviewer', required=True, help='Who decided')
    parser.add_argument(
        UNRESOLVE_ARG,
        action='store_true',
        help='Reopen a manually resolved flag instead, so it returns to the report',
    )
    parser.add_argument('--yes', dest='assume_yes', action='store_true', help='Skip the confirmation prompt')
    args = parser.parse_args()
    return main(
        dataset=args.dataset,
        sg_ids=args.sg_ids,
        category=args.category,
        reason=args.reason,
        reviewer=args.reviewer,
        unresolve=args.unresolve,
        assume_yes=args.assume_yes,
    )


if __name__ == '__main__':
    raise SystemExit(cli_main())
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest test -q`
Expected: PASS, whole suite.

- [ ] **Step 5: Commit**

```bash
git add src/rd_qc/scripts/resolve_somalier_flag.py test/test_resolve_somalier_flag.py
git commit -m "feat: add the resolve_somalier_flag CLI"
```

---

### Task 7: Entry point and documentation

**Files:**
- Modify: `pyproject.toml:43-44`
- Modify: `README.md` (the `## Flags` section, and the directory structure listing)

- [ ] **Step 1: Add the console script**

In `pyproject.toml`, under `[project.scripts]`:

```toml
[project.scripts]
run_workflow = 'rd_qc.run_workflow:cli_main'
resolve_somalier_flag = 'rd_qc.scripts.resolve_somalier_flag:cli_main'
```

- [ ] **Step 2: Check it resolves**

Run: `uv run resolve_somalier_flag --help`
Expected: the argparse help, listing `--dataset`, `--sg-ids`, `--category` with its three choices, `--reason`, `--reviewer`, `--unresolve`, `--yes`.

- [ ] **Step 3: Document it in the README**

Append to the `## Flags` section, after the paragraph ending "so the report can show a backlog of past findings alongside current ones":

````markdown
A flag that is real but accepted, such as a pedigree known to be wrong that will not be corrected, can be resolved by hand:

```bash
resolve_somalier_flag --dataset my-dataset --sg-ids CPG001 CPG002 \
    --category relatedness_mismatch --reason "pedigree known wrong" --reviewer ef
```

Take the sequencing group IDs off the report row in either order. A manually resolved flag stays resolved even while the checks keep measuring the finding, and it is left out of the report entirely rather than shown as resolved history. The resolution is bound to that exact finding, so if the genotypes later say something different about the same pair, the new finding is reported as usual. `--unresolve` reverses it.
````

Add `resolve_somalier_flag.py` and `flag_store.py` to the directory structure listing in the same file:

```
src
└── rd_qc
    ├── run_workflow.py
    ├── config_template.toml
    ├── flag_store.py
    ├── stages.py
    ├── utils.py
    ├── jobs
    │   ├── generate_somalier.py
    │   ├── relate.py
    │   └── somalier_flags_report.py
    ├── scripts
    │   ├── check_pedigree.py
    │   ├── check_self_relatedness.py
    │   ├── record_somalier_flags.py
    │   ├── resolve_somalier_flag.py
    │   └── somalier_flags_report.py
    └── templates
        └── somalier_flags_overview.html.jinja
```

And add a line to the `## Key Components` section, after the `utils.py` paragraph:

```markdown
`flag_store.py` is the only module that reads or writes the `somalier_flags` list on a sequencing group's meta. The mutation replaces the whole list, so both the pipeline's reconciler and `resolve_somalier_flag` go through it.
```

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml README.md
git commit -m "docs: document the manual flag resolution CLI"
```

---

### Task 8: Verification and version bump

**Files:**
- Modify: `pyproject.toml`, `Dockerfile`, `README.md` (all three via bump-my-version)

- [ ] **Step 1: Run the whole suite**

Run: `uv run pytest test -q`
Expected: PASS, no skips.

- [ ] **Step 2: Check coverage of the new code**

Run: `uv run pytest test --cov=src/rd_qc --cov-report=term-missing -q`
Expected: `flag_store.py`, `resolve_somalier_flag.py` and the changed branches in `record_somalier_flags.py` at 80% or better. The only expected misses are `cli_main` and the `input()` call in `confirmed`, both of which are argument plumbing.

- [ ] **Step 3: Lint and format**

`ruff` is not a project dependency, so plain `uv run ruff` fails with "Failed to spawn". The version is pinned to match `.pre-commit-config.yaml`, because a newer ruff turns on `PLR0917` and reports 11 pre-existing findings in files this work does not touch.

Run: `uv run --with ruff==0.15.19 ruff check src test && uv run --with ruff==0.15.19 ruff format --check src test`
Expected: `All checks passed!` and `N files already formatted`. Fix anything reported with `uv run --with ruff==0.15.19 ruff format src test` and re-run.

- [ ] **Step 4: Bump the version**

A new user-facing capability, so minor: 0.2.4 to 0.3.0.

Run: `uv run bump-my-version bump minor`
Expected: `pyproject.toml`, `Dockerfile` and `README.md` updated and committed as `Bump version: 0.2.4 → 0.3.0`.

- [ ] **Step 5: Review the whole diff**

Run: `git diff main...HEAD`

Check specifically:
- `record_somalier_flags.py` has no leftover `query` import or GraphQL literal
- the manual branch sits **above** the `compare_*` branch in all three reconcilers, not below it
- no `manually_resolved` handling leaked into `check_pedigree.py` or `check_self_relatedness.py`, which must keep measuring and emitting everything

---

## Out of scope

Bulk resolution, a `list` subcommand for finding flags outside the report, and any surfacing of reasons or reviewers in the report or Slack. See the spec's own out-of-scope section.
