# Qodex.summary

## Task
Migrate all remaining V-LiSEMOD database-backed routes/pages to RANDY strict backend mode.

## Original Goal
The user wants every page, not just the homepage workflow, to use RANDY/background apps for database access so Heroku does not fail on missing local SQLite tables.

## Assumptions
- Heroku V-LiSEMOD should run with `VLISMOD_DATA_BACKEND=randy`.
- The production RANDY base is `https://randy.rove-vernier.ts.net/backup/vlismod`.
- The production database is hosted on RANDY at `/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`.
- Heroku may still generate request-time artifacts locally, including ligand SVGs, charts, cached SDFs, and `.pml` files.
- Strict `randy` mode must not silently fall back to local SQLite.
- Debug-only routes may remain local if they are not used by the production UI flow.

## Files Inspected
- `app.py`
- `RANDY/app.py`
- `RANDY/vlismod_data_routes.py`
- `templates/compare_ligands.html`
- `templates/query_protein_virus.html`
- `templates/protacability_assessment.html`
- `templates/display_images.html`
- `templates/ligand_query.html`
- `README.md`
- `docs/DEPLOYMENT.md`
- Remote runtime file: `/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`

## Files Changed
- `app.py`
- `RANDY/vlismod_data_routes.py`
- `README.md`
- `docs/DEPLOYMENT.md`
- `Qodex.summary.md`

## Files Created
- `tools/check_vlismod_randy_routes.py`

## Implementation Summary
Completed an app-wide route audit and migrated the remaining production route groups that still performed local SQLite reads in user-facing flows.

New RANDY-backed route families now cover:
- Compare Ligands data loading
- ligand synonym/info/option lookup
- interaction comparison and chart data loading
- Protein Query filter and export data loading
- PROTACability filter/search/detail/export data loading
- ligand image generation input data
- PyMOL session input data

The main app now uses RANDY for these routes in strict `VLISMOD_DATA_BACKEND=randy` mode and preserves local fallback only for `local` or `auto`.

Also removed a hidden local-DB dependency from the ligand-instance mapping flow used by ligand SDF/viewer routes. In strict `randy` mode, that flow no longer opens local SQLite just to enrich diagnostics.

Added `RANDY_API_TIMEOUT_SECONDS` support so larger RANDY-backed payloads are not forced through a short hardcoded client timeout.

## Route Audit Summary
- Static/no DB pages:
  - landing/info pages
  - template-only navigation pages
  - coordinate-file and cache-serving routes that use local files, not SQLite
- Already RANDY-backed before this audit:
  - `/get_viruses`
  - `/get_pdb_codes/<virus_name>`
  - `/get_ligands/<pdb_code>`
  - `/check_functional_groups/<pdb_code>`
  - `/get_ligands_list`
  - `/get_viruses_by_ligand/<ligand_code>`
  - `/get_pdb_residue_by_ligand/<ligand_code>`
  - `/get_pdb_mapping/<ligand_code>`
  - `/get_sasa_chains/<pdb_code>/<ligand_name>`
  - `/generate_ligand_images`
  - `/generate_pymol_session`
- Newly migrated in this audit:
  - `/get_ligands_with_synonyms`
  - `/get_ligand_info/<ligand_code>`
  - `/get_ligand_options`
  - `/get_smiles_svg/<ligand_id>`
  - `/get_atom_count/<ligand_id>`
  - `/highlight_atoms`
  - `/generate_charts`
  - `/compare_ligand_interactions`
  - `/get_virus_names_list_distinct`
  - `/get_protein_types_list_distinct`
  - `/get_pdbs_for_virus_protein`
  - `/export_data_to_excel`
  - `/query_protein_virus_page`
  - `/protacability_page`
  - `/api/protacability/filter_options`
  - `/api/protacability/search`
  - `/api/protacability/detail/<pdb_code>/<chain_id>`
  - `/api/protacability/structure_detail/<pdb_code>`
  - `/api/protacability/protein_detail`
  - `/api/protacability/target_detail`
  - `/api/protacability/export`
  - `/api/ligand_instance_sdf_url/<pdb_code>/<ligand_code>`
  - `/api/ligand_instance_sdf/<pdb_code>/<ligand_code>.sdf`
- Deferred or intentionally local:
  - `/api/debug/*` routes
  - routes that inspect local coordinate assets for diagnostics
  - helper code paths used only for debug payloads

## Key Decisions
- Heroku must not read production SQLite directly because the production database lives on RANDY, not on the dyno filesystem.
- RANDY owns database reads; V-LiSEMOD owns request-time rendering and file generation.
- Strict `randy` mode should fail loudly with remote errors instead of silently creating or reading an empty local SQLite file.
- Local SVG, chart, coordinate-cache, and PyMOL artifact generation can remain on Heroku because they are derived outputs, not authoritative data storage.
- Debug-only routes were left local where reasonable because they are not part of the production user flow and do not need immediate RANDY migration.

## Commands Run
- `python -m py_compile app.py RANDY/app.py RANDY/vlismod_data_routes.py tools/check_vlismod_randy_routes.py`
- `rg`/`grep` route and SQLite audit commands over `app.py`, `templates/`, and `static/`
- local Flask test-client checks in `local` mode
- local Flask test-client checks in strict `randy` mode
- remote sync:
  - `scp RANDY/vlismod_data_routes.py randy:/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`
- remote Gunicorn reload:
  - `ssh randy 'kill -HUP <backup_receiver master pid>'`
- remote Flask test-client checks against `backup_receiver.app:APP`
- remote schema inspection with remote Python `sqlite3`
- direct HTTPS `curl`/`requests` checks against `https://randy.rove-vernier.ts.net/backup/vlismod`

## Validation Results
- Syntax:
  - `python -m py_compile ...` passed.

- Local `local`-mode regression checks passed:
  - `/get_ligand_info/DR7` -> `200`
  - `/compare_ligand_interactions` -> `200`
  - `/api/protacability/filter_options` -> `200`
  - `/api/protacability/search?view=chain&page=1&page_size=5` -> `200`
  - `/api/protacability/detail/1A43/A` -> `200`
  - `/api/ligand_instance_sdf_url/3EKY/DR7?auth_chain=A&auth_seq_id=100` -> `200`

- RANDY blueprint validation in local workspace passed for the newly added endpoints, including:
  - `/backup/vlismod/ligands/with-synonyms`
  - `/backup/vlismod/ligand-info`
  - `/backup/vlismod/ligand-options`
  - `/backup/vlismod/interaction-records`
  - `/backup/vlismod/ligand-interactions/compare`
  - `/backup/vlismod/virus-proteins/virus-names`
  - `/backup/vlismod/virus-proteins/protein-types`
  - `/backup/vlismod/virus-proteins/pdbs`
  - `/backup/vlismod/export-data`
  - `/backup/vlismod/protacability/source`
  - `/backup/vlismod/protacability/raw-table`

- Remote live RANDY runtime validation on the RANDY host passed via Flask test client against `backup_receiver.app:APP`:
  - `/backup/vlismod/ligands/with-synonyms` -> `200`
  - `/backup/vlismod/ligand-info?ligand_code=DR7` -> `200`
  - `/backup/vlismod/ligand-options` -> `200`
  - `/backup/vlismod/virus-proteins/virus-names` -> `200`
  - `/backup/vlismod/virus-proteins/protein-types` -> `200`
  - `/backup/vlismod/protacability/source` -> `200`

- Strict `randy` mode guardrail behavior was confirmed earlier in the app:
  - when RANDY routes are missing or remote calls fail, the app returns visible remote errors instead of touching local SQLite.

- Important limitation:
  - direct external HTTPS validation of some newly added heavy RANDY endpoints from this Codex environment remained inconsistent and often timed out, even after the live route file was synced and internal RANDY-host validation passed.
  - This means the endpoint logic is proven on the running RANDY Flask app itself, but full external-path validation for every new route still needs one more post-deploy smoke test from Heroku or another network vantage point.

## Known Issues
- External/public HTTPS access to the newer heavy RANDY endpoints was not consistently verifiable from this environment; internal live RANDY validation passed, but public-path validation remains incomplete for some routes.
- Debug routes under `/api/debug/*` still use local SQLite and were intentionally not migrated.
- The ligand-instance mapping routes no longer require local SQLite in strict `randy` mode, but they still depend on local coordinate assets or RCSB CIF/SDF fetches for actual structural resolution.
- Any routes not exercised by the UI or not included in the major audited route groups should still get one production smoke pass after deploy.

## Manual Verification
1. Redeploy V-LiSEMOD with:
   - `VLISMOD_DATA_BACKEND=randy`
   - `VLISMOD_BACKUP_URL=https://randy.rove-vernier.ts.net/backup/vlismod`
   - `RANDY_API_TOKEN=<secret>`
   - `RANDY_API_TIMEOUT_SECONDS=45`
2. Run direct RANDY checks:
   - `curl -i -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/db-health"`
   - `curl -i -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/ligands/with-synonyms"`
   - `curl -i -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/virus-proteins/virus-names"`
   - `curl -i -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/protacability/source"`
3. Run the route checker:
   - `python tools/check_vlismod_randy_routes.py`
4. Browser smoke test:
   - homepage workflow
   - Compare Ligands page
   - Ligand Indexer page
   - Protein Query page
   - PROTACability page
   - ligand image generation
   - PyMOL session generation
5. Watch Heroku logs for any remaining `sqlite3.OperationalError: no such table` messages.

## Suggested Next Prompt
If the post-deploy smoke test shows any remaining gaps, use:

“Audit the remaining debug/structure-helper routes and any coordinate-asset-dependent viewer endpoints, then migrate or explicitly production-disable any path that still touches local SQLite in strict `VLISMOD_DATA_BACKEND=randy` mode.”
