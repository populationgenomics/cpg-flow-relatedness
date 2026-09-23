# CPG-Flow - Rare Disease QC

Version 0.3.0

A CPG workflow for sample identity and relatedness QC using [Somalier](https://github.com/brentp/somalier), built on the [cpg-flow](https://github.com/populationgenomics/cpg-flow) pipeline framework.

## Purpose

This workflow verifies sample identity and pedigree consistency for rare disease cohorts. It uses Somalier to extract genomic fingerprints from sequencing data (CRAMs or gVCFs), then runs cross-sample identity checks and pedigree validation. Findings are recorded in Metamist as structured flags, reported to Slack, and published as an interactive HTML report.

It is intended to catch sample swaps, contamination, and pedigree errors early in the analysis pipeline, before downstream variant interpretation in [seqr](https://github.com/populationgenomics/seqr).

## Workflow Overview

1. Query Metamist for all sequencing groups in the project and their existing Somalier fingerprints
2. Generate Somalier fingerprints for any sequencing groups that are missing them, selecting the best available source file (CRAM > gVCF, configurable)
3. Run self-relatedness checks for participants with multiple sequencing groups (e.g. short-read + long-read), comparing fingerprints to verify they belong to the same individual
4. Run a pedigree check across all sequencing groups, comparing the pedigree Metamist holds against what the genotypes measure
5. Record the resulting flags against each sequencing group's `meta` in Metamist, resolving any that no longer apply
6. Render every recorded flag into a single report, grouped by family

## Flags

Three categories of flag are recorded, each keyed by its own attributes:

- `sex_inference_mismatch` — the inferred sex disagrees with the sex on file
- `self_relatedness_mismatch` — two sequencing groups from one participant are not related enough to be the same person
- `relatedness_mismatch` — the measured relatedness disagrees with the pedigree

Flags persist across runs. A flag that reappears is updated in place, and one that no longer applies is marked resolved with a date rather than deleted, so the report can show a backlog of past findings alongside current ones. A run that changes nothing about a sequencing group's flags writes nothing, so the meta history only records real changes.

A flag that is real but accepted, such as a pedigree known to be wrong that will not be corrected, can be resolved by hand:

```bash
resolve_somalier_flag --dataset my-dataset --sg-ids CPG001 CPG002 \
    --category relatedness_mismatch --reason "pedigree known wrong" --reviewer ef
```

Take the sequencing group IDs off the report row in either order. A manually resolved flag stays resolved even while the checks keep measuring the finding, and it is left out of the report and the Slack summary counts entirely, rather than shown as resolved history. The record is closed at the point the curator resolved it: later runs leave its measurements alone rather than refreshing them. The resolution is bound to that exact finding, so if the genotypes later say something different about the same pair, the new finding is reported as usual. Run it locally, against your own Metamist credentials. If more than one flag matches the category and sequencing group IDs it refuses to guess, listing the candidates instead. `--unresolve` reverses a resolution.

## Relatedness inferences

Relationships are inferred from what Somalier measured, rather than from its `somalier relate --infer` pedigree reconstruction. Each pair's kinship coefficient places it in a degree band (identical, parent-child, siblings, second-, third-degree, unrelated), with IBS0 separating parent-child from full siblings. That measured degree is then compared against the degrees the recorded pedigree allows.

`--infer` is not used because it assumes high quality pairs with both parents present. Pedigrees that record only one parent send it well outside that envelope, where it rewrites family IDs and invents parent links that a pedigree check then reads as truth.

Not every disagreement is treated as an error:

- **Conflicts** are cases where the pedigree and the genotypes genuinely contradict each other. These lead the report.
- **Refinements** are cases where the pedigree is merely less specific than the genotypes, most often two people recorded in one family with no path between them because a parent is missing from the database. Each is a candidate pedigree correction, and they usually outnumber conflicts, so they get their own de-emphasised section.

Two cases are deliberately left unflagged. Third-degree relatedness between families is not reported at all, because at that distance the measurement is indistinguishable from cohort background. And related parents are not reported when a child of their union is marked as consanguineous in Metamist, since the pedigree already accounts for the relatedness.

## Directory Structure

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

## Key Components

`config_template.toml` contains the default configuration for the workflow, including source file priority for fingerprint extraction, Somalier reference sites, and Slack notification settings.

`stages.py` defines four `DatasetStage` subclasses that run once per dataset in the cohort:

- `GenerateMissingSomalierFingerprints` — extracts fingerprints only for sequencing groups that don't already have them
- `SomalierSelfCheck` — compares fingerprints across sequencing groups for participants with two or more SGs
- `SomalierPedigreeCheck` — compares measured relatedness against the expected pedigree
- `GenerateSomalierFlagsReport` — renders every flag recorded in Metamist into the HTML report

`utils.py` contains the flag dataclasses, the relatedness inference and verdict logic, and utility functions for querying Metamist and building PED files.

`flag_store.py` is the only module that writes the `somalier_flags` list on a sequencing group's meta. The mutation replaces the whole list rather than patching entries, so a second copy of that write is a second chance to drop every flag on a sequencing group, and both the pipeline's reconciler and `resolve_somalier_flag` go through this one. Reading is less dangerous, and the reconciler and the report still read the meta key directly.

`jobs/` contains the Hail Batch job builders. `scripts/` contains the post-processing scripts those jobs run.

## Running the tests

```bash
uv run pytest test/ -q
```

The suite runs entirely offline against synthetic Somalier outputs and mock fixtures, so it needs no database or network access.

## Rendering the report locally

Two scripts under `testing_scripts/` render the report without running the workflow.

For layout and styling work, render the mock fixtures. This needs no credentials:

```bash
uv run python testing_scripts/render_mock_somalier_report.py
```

To check behaviour against a real dataset, download `samples.tsv` and `pairs.tsv` from a `SomalierPedigreeCheck` run into `local_data/<dataset>/`, then:

```bash
SM_ENVIRONMENT=production uv run python testing_scripts/local_pedigree_report.py \
    --input-dir local_data/<dataset> --dataset <dataset>
```

This rebuilds the expected pedigree from Metamist, re-derives the flags from the TSVs and writes `report.html`. It is read-only: no SG meta is modified, no analysis is registered and nothing is sent to Slack, so it is safe to point at a production dataset. The Metamist responses are cached alongside the inputs, so later renders can add `--offline` and skip the network. The first run for a dataset has to be online, since the expected pedigree comes from Metamist.
