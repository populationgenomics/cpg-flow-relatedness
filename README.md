# cpg-flow-pipeline-template

A template repository to use as a base for CPG workflows using the [cpg-flow](https://github.com/populationgenomics/cpg-flow) pipeline framework.

Current version: 0.1.1

> **Upstream**: This template builds on the conventions from [cpg-python-template-repo](https://github.com/populationgenomics/cpg-python-template-repo), adding cpg-flow-specific architecture.

## Quick Start

1. Create a new repo from this template
2. Rename `workflow_name` and `cpg_flow_<tool>` everywhere — directory name under `src/`, `pyproject.toml`, `Dockerfile`, imports, `get_version.py`
3. Update `pyproject.toml`: set name, description, add tool-specific dependencies, and check version here and in the Dockerfile and README.md (default is `0.1.0`)
4. **Update `Dockerfile`**: Set `WORKDIR`, `ENV VERSION`
5. Update `config_template.toml`: replace placeholder `[tool_name]` config with your tool's settings
6. Update `run_workflow.py`: set workflow name and wire up your top-level stage(s)
7. ** Write your stages, jobs, and scripts **
8. **Update `.pre-commit-config.yaml`**: Ensure ruff and mypy point to the correct package name
9. Set up GitHub secrets for CI/CD
10. Verify: `pip install .[test]`, `pre-commit run --all-files`, `pytest test`, `docker build .`
11. **Verify quotes**: Ensure bump-my-version and ruff agree on quote style (`'` single quotes — ruff enforces this via `Q000` ignore + `quote-style = 'single'`)

## Directory Structure

```
src
├── workflow_name
│   ├── __init__.py
│   ├── config_template.toml
│   ├── jobs
│   │   └── LogicForAStage.py
│   ├── run_workflow.py
│   ├── stages.py
│   └── utils.py
```

`workflow_name` occurs in a number of places ([pyproject.toml](pyproject.toml), [src](src), and the workflow name in the template
config file). Crucially it also appears in [the image builder workflow](.github/workflows/get_version.py), which
determines the names for images built from this repository. It is intended that you remove this generic placeholder
name, and replace it with the name of your workflow.

`stages.py` contains Stages in the workflow, with the actual logic imported from files in `jobs`.

`stages.py` also links to the Pipeline Naming Conventions document, containing a number of recommendations for naming
Stages and other elements of the workflow.

`config_template.toml` is a template, indicating the settings which are mandatory for the pipeline to run. In
production-pipelines, many of these settings were satisfied by the cpg-workflows or per-workflow default TOML files. If
a pipeline is being migrated from production-pipelines, the previous default config TOML would be a better substitute.
