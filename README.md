# CPG-Flow - Rare Disease QC

Version 0.1.6

A CPG workflow for sample identity and relatedness QC using [Somalier](https://github.com/brentp/somalier), built on the [cpg-flow](https://github.com/populationgenomics/cpg-flow) pipeline framework.

## Purpose

This workflow verifies sample identity and pedigree consistency for rare disease cohorts. It uses Somalier to extract genomic fingerprints from sequencing data (CRAMs, gVCFs, or VCFs), then runs cross-sample identity checks and pedigree validation. Results are reported to Slack and published as interactive HTML reports.

It is intended to catch sample swaps, contamination, and pedigree errors early in the analysis pipeline, before downstream variant interpretation in [seqr](https://github.com/populationgenomics/seqr).

## Workflow Overview

1. Query Metamist for all sequencing groups in the project and their existing Somalier fingerprints
2. Generate Somalier fingerprints for any sequencing groups that are missing them, selecting the best available source file (VCF > gVCF > CRAM, configurable)
3. Run cross-type identity checks for participants with multiple sequencing groups (e.g. short-read + long-read), comparing fingerprints to verify they belong to the same individual
4. Run a pedigree check across all sequencing groups using the project pedigree from Metamist, validating expected relatedness between family members

## Directory Structure

```
src
└── rd_qc
    ├── __init__.py
    ├── run_workflow.py
    ├── config_template.toml
    ├── stages.py
    ├── utils.py
    ├── jobs
    │   ├── generate_somalier.py
    │   └── relate.py
    └── scripts
        ├── __init__.py
        ├── check_pedigree.py
        └── check_self_relatedness.py
```

## Key Components

`config_template.toml` contains the default configuration for the workflow, including source file priority for fingerprint extraction, Somalier reference sites, and Slack notification settings.

`stages.py` defines three `DatasetStage` subclasses that run once per dataset in the cohort:
- `GenerateMissingSomalierFingerprints` — extracts fingerprints only for sequencing groups that don't already have them
- `RunCrossTypeIdentityChecks` — compares fingerprints across sequencing groups for participants with two or more SGs
- `SomalierPedigreeCheck` — validates relatedness against the expected pedigree

`utils.py` contains utility functions for querying Metamist (sequencing groups, analyses, pedigrees), managing the `SomalierIndex` data structure, selecting source files for extraction, and building PED files.

`jobs/` contains the Hail Batch job builders for Somalier extract and relate commands.

`scripts/` contains post-processing scripts that run as Hail Batch jobs, including Slack reporting for identity and pedigree check results.
