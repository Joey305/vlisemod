# Maintenance Guide

## Generated Folders That May Exist

During normal use or enrichment workflows, the following generated folders may appear or be populated:

- `output_csvs/`
- `pml_sessions/`
- `temp/`
- `users_info/`
- `static/charts/`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- `static/ligand_images/`
- `PDB_FILES/`

Treat these as runtime or pipeline artifacts unless a deliberately curated sample is being prepared.

## Cache Cleanup

Safe maintenance targets usually include:

- `static/charts/`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- temporary export folders such as `output_csvs/` and `temp/`

Before cleanup, confirm nothing still depends on the generated files for an active demo or local run.

## Export / Output Folders

- `output_csvs/` is used for dataset export assembly.
- `pml_sessions/` stores generated PyMOL-related outputs.
- temporary download/export bundles may be created and removed during request workflows.

If exports start failing, check for path existence, permissions, and leftover partial files.

## Static Generated Assets

Likely generated static assets include:

- chart images
- ligand images
- cached coordinates for viewers
- cached ligand SDF files

These should usually remain ignored by Git and can often be regenerated.

## Local Database Handling

- `viral_data.db` is the main app database and should be handled as a provisioned local/runtime asset.
- `users.db` is deprecated and no longer used by the application.
- before replacing or regenerating a local database, make a separate backup outside routine Git operations.

## Files / Folders That Should Stay Gitignored

Confirmed ignore targets already cover the main risk areas:

- `.env`
- `users.db`
- `viral_data.db`
- `PDB_FILES/`
- `output_csvs/`
- `pml_sessions/`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- `static/charts/`
- `static/ligand_images/`
- model/checkpoint folders
- `__pycache__/`
- `.DS_Store`

Keep these ignored unless a very small redacted sample is intentionally prepared.

## Updating Documentation After Route / Page Changes

Whenever public routes, page labels, or navigation behavior changes:

1. update `README.md` only if the public-facing product story changes,
2. update `docs/APP_GUIDE.md` for page/module changes,
3. update `docs/DEVELOPER_NOTES.md` for route or feature-flag changes,
4. update `docs/DEPLOYMENT.md` if env vars or runtime behavior changes.

## Recommended Pre-Push Checklist

1. Confirm no local databases, caches, exports, or model weights are staged.
2. Confirm `.env` and any secret-bearing files remain ignored.
3. Confirm README and docs links work with relative GitHub paths.
4. Confirm public wording does not overclaim PROTACability results.
5. Confirm any route references still exist in the app.
6. Confirm optional-module notes still match the current feature flags.

## Troubleshooting Common App Issues

### App does not start

- verify dependencies from `requirements.txt` are installed
- watch for missing packages such as `seaborn`, `rdkit`, or `biopython`

### Query pages load but data is empty

- verify `viral_data.db` exists and contains the expected tables
- verify the selected filters actually map to current database contents

### PROTACability page shows unavailable data

- import the CSV-backed PROTACability tables using `TOOLS/import_protacability_data.py`
- verify the required CSV files exist in `PDB_FILES/`

### Drug GPT behaves inconsistently

- verify `SHOW_DRUG_GPT_NAV`, `ENABLE_DRUG_GPT`, and `ENABLE_LOCAL_LLM`
- remember that `/drugapp/` intentionally returns a disabled-state response when the feature is off

### Generated outputs fail

- verify write permissions for generated folders
- verify the source tables or structure files required by that workflow exist
