# Qodex.summary

## Task
Set up authenticated RANDY API backend support for V-Lismod database access.

## Original Goal
Offload V-LiSEMOD database-backed lookup operations to RANDY so V-LiSEMOD does not require a local `viral_data.db` backup in production, while keeping the app deployable with Heroku-style environment-variable configuration and preserving local development fallback behavior.

## Assumptions
- RANDY and V-LiSEMOD are separate Flask apps that may be deployed independently.
- Local testing uses `/Users/jxs794/Documents/VLISEMOD/viral_data.db` unless `VLISMOD_DB_PATH` is set differently.
- The first migration batch should cover only the simple JSON lookup routes already used by templates and frontend JavaScript.
- Complex chart generation, ligand rendering, PyMOL/session generation, interaction comparisons, SASA highlighting, and PROTACability routes remain local for now.

## Files Inspected
- `app.py`: identified SQLite-backed routes and preserved response shapes.
- `RANDY/app.py`: checked Flask app creation and blueprint registration style.
- `RANDY/e3_data_routes.py`: reused the existing bearer-token and JSON error-handling pattern.
- `requirements.txt`: confirmed `requests` was already available.
- `README.md`: updated high-level backend configuration notes.
- `docs/DEPLOYMENT.md`: updated deployment and Heroku-style environment-variable guidance.
- `Dockerfile`: confirmed no current Heroku runtime wiring lived there.
- `setup.sh`: confirmed it was unrelated to the API migration.

## Files Changed
- `app.py`: added RANDY client helpers and switched the supported lookup routes to local/RANDY/auto backend dispatch.
- `RANDY/app.py`: registered the new V-LiSEMOD blueprint and made imports work in both package-style and local script execution.
- `README.md`: added a short note about remote RANDY-backed data mode.
- `docs/DEPLOYMENT.md`: documented RANDY/V-LiSEMOD environment variables and local verification flow.

## Files Created
- `RANDY/vlismod_data_routes.py`: authenticated RANDY blueprint for V-LiSEMOD data access.
- `Qodex.summary.md`: this summary document.

## Implementation Summary
Created a new RANDY blueprint at `/api/vlismod` that reads the database path from `VLISMOD_DB_PATH`, protects routes with `Authorization: Bearer <token>` using `VLISMOD_API_TOKEN`, and returns JSON-only responses with structured error handling. Added health and database-health endpoints plus the first batch of V-LiSEMOD lookup endpoints.

Updated the main V-LiSEMOD Flask app to support `VLISMOD_DATA_BACKEND=local|randy|auto`, using `RANDY_API_BASE_URL` and `RANDY_API_TOKEN` for authenticated HTTP calls to RANDY. In `auto` mode, supported routes fall back to the local SQLite queries if the RANDY request fails. Public V-LiSEMOD route URLs and response shapes were preserved.

## Key Decisions
- Token auth is environment-variable based so secrets stay out of source control and can be configured cleanly on Heroku.
- Local SQLite fallback was preserved to avoid breaking current development and to keep rollout risk low.
- Only the simple lookup routes were migrated first because they have stable JSON contracts and low coupling to local file generation or heavy processing.

## Commands Run
- Searched for SQLite usage, `viral_data.db` references, blueprint registration, and authorization handling with `rg`.
- Inspected relevant file sections with `sed`.
- Ran `python -m py_compile app.py RANDY/app.py RANDY/e3_data_routes.py RANDY/vlismod_data_routes.py`.
- Ran import checks for both Flask apps with short Python snippets.
- Ran Flask `test_client()` smoke checks for RANDY auth and V-LiSEMOD local/proxied lookup routes.
- Ran a targeted secret-safety search over changed files.

## Validation Results
- Syntax compilation passed for `app.py`, `RANDY/app.py`, `RANDY/e3_data_routes.py`, and `RANDY/vlismod_data_routes.py`.
- Main V-LiSEMOD app import check passed.
- RANDY app import check passed when loaded with the `RANDY/` directory on `sys.path`, matching normal script execution.
- RANDY `test_client()` checks passed for:
  - `/api/vlismod/health`
  - unauthorized `/api/vlismod/viruses` returning `401`
  - authenticated `/api/vlismod/db-health`
  - authenticated `/api/vlismod/viruses`
- Main V-LiSEMOD `test_client()` checks passed in local mode for:
  - `/get_viruses`
  - `/get_ligands_list`
  - `/check_functional_groups/6Y2F`
- Main V-LiSEMOD proxy-mode checks passed with a mocked RANDY HTTP layer for:
  - `/get_viruses`
  - `/get_ligands_list`
  - `/get_pdb_mapping/NAG`

## Known Issues
- Complex routes that build charts, images, PyMOL sessions, interaction payloads, or PROTACability outputs still query local SQLite directly.
- There is still no root `Procfile`; deployment is documented but not fully packaged for Heroku within this repository.
- The placeholder token string used in docs for manual verification should be replaced with a real secret only in environment configuration, never in code.

## Manual Verification
1. Start RANDY locally with `VLISMOD_API_TOKEN` and `VLISMOD_DB_PATH`.
2. Call `/api/vlismod/db-health` with `Authorization: Bearer <token>`.
3. Start V-LiSEMOD with `VLISMOD_DATA_BACKEND=randy` or `auto`, plus `RANDY_API_BASE_URL` and `RANDY_API_TOKEN`.
4. Open pages that use virus, PDB, ligand, mapping, or SASA lookup routes.
5. Confirm those pages load without V-LiSEMOD directly needing a local database backup in strict RANDY mode.

## Suggested Next Prompt
Migrate the next batch of V-LiSEMOD database-backed routes to RANDY by adding authenticated API endpoints and proxy support for ligand interaction payloads, ligand comparison data, SASA highlighting support queries, `get_ligand_options`, `get_ligands_with_synonyms`, `get_ligand_info`, and the first safe subset of PROTACability JSON endpoints, while keeping chart/image generation local unless the data endpoint is straightforward to expose.
