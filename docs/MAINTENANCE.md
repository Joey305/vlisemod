# Maintenance Guide

## Generated Folders That May Exist

During normal use or enrichment workflows, the following generated folders may appear or be populated:

- `output_csvs/`
- `pml_sessions/`
- `temp/`
- `tmp/`
- `users_info/`
- `static/charts/`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- `static/ligand_images/`
- `PDB_FILES/`

Treat these as runtime or pipeline artifacts unless a deliberately curated public sample is being prepared.

## Public GitHub Hygiene

Before making the repository public or sharing a branch:

1. confirm local databases, caches, exports, logs, and archives are not staged,
2. confirm `.env` and secret-bearing local config remain ignored,
3. confirm model weights and checkpoints remain outside version control,
4. confirm documentation does not expose private infrastructure details.

If risky files are already tracked, remove them from Git history or at least untrack them in a follow-up cleanup step before public release.

## Cache and Output Cleanup

Safe cleanup targets usually include:

- `static/charts/`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- temporary export folders such as `output_csvs/`, `temp/`, and `tmp/`

Before cleanup, confirm no active demo or local run still depends on the generated files.

## Local Database Handling

- `viral_data.db` is a provisioned runtime asset, not source code.
- `users.db` is deprecated and should not be part of the public workflow.
- before replacing or regenerating a local database, keep a separate backup outside normal Git operations.

## Files and Folders That Should Stay Ignored

High-priority ignore targets include:

- `.env`
- `viral_data.db`
- `users.db`
- `PDB_FILES/`
- `output_csvs/`
- `pml_sessions/`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- `static/charts/`
- `static/ligand_images/`
- `models/`
- `checkpoints/`
- `__pycache__/`
- `.DS_Store`
- `*.log`
- large archives such as `*.zip`, `*.tar`, and `*.tar.gz`

## Updating Docs After App Changes

Whenever public routes, page labels, feature flags, or backend behavior change:

1. update `README.md` if the public-facing story changed,
2. update `docs/APP_GUIDE.md` for module and workflow changes,
3. update `docs/DEPLOYMENT.md` for env vars or runtime changes,
4. update public manuscript or citation materials once they are added to the documentation tree.

## Troubleshooting Common App Issues

### App does not start

- verify dependencies from `requirements.txt` or `environment.yml`
- watch for missing packages such as `rdkit`, `matplotlib`, `seaborn`, or `biopython`

### Query pages load but data is empty

- verify `viral_data.db` exists and contains the expected tables
- verify filters map to current database contents

### PROTACability page shows unavailable data

- import the PROTACability tables using `TOOLS/import_protacability_data.py`
- verify the required CSV inputs exist in `PDB_FILES/`

### Optional assistant behaves unexpectedly

- verify `SHOW_DRUG_GPT_NAV`, `ENABLE_DRUG_GPT`, and `ENABLE_LOCAL_LLM`
- remember that `/drugapp/` can intentionally return a disabled-state response when the feature is off

### Generated outputs fail

- verify write permissions for generated folders
- verify the source tables or structure files required by that workflow exist
