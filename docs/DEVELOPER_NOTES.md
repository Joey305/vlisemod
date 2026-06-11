# Developer Notes

## Flask App Structure

The main application lives in `app.py` and exposes most public routes directly from a single Flask app instance. An optional Drug GPT assistant module lives in `DRUGapp.py` and is registered as the `dp` blueprint at `/drugapp` only when `ENABLE_DRUG_GPT` is enabled.

The current public workflow is no-login by default. Authentication routes, `users.db`, and Flask-Login are not part of the active runtime path.

## High-Level Route Areas

Key page routes currently include:

- `/`
- `/about`
- `/ligand_indexer`
- `/compare_ligands`
- `/query_protein_virus_page`
- `/protacability_page`
- `/healthz`

Key API areas currently include:

- ligand lookup and mapping endpoints
- image and PyMOL-session generation endpoints
- export endpoints
- PROTACability filter, search, detail, and export endpoints

## Feature Flags and Optional Modules

Current flags in `app.py`:

- `SHOW_DRUG_GPT_NAV`
- `ENABLE_DRUG_GPT`
- `ENABLE_LOCAL_LLM`

If Drug GPT is disabled, `app.py` still serves disabled-state routes for `/drugapp/` and `/drugapp/query` so the UI can fail gracefully instead of hard-404ing.

## RANDY and Local Data Paths

The app supports:

- direct local SQLite access through `viral_data.db`, and
- RANDY-backed access for selected route groups when `VLISMOD_DATA_BACKEND` is set to `randy` or `auto`.

This is important for both documentation and debugging because some user-facing workflows may look identical while drawing data from different backends.

## Template and Navigation Conventions

- shared shell: `templates/base.html`
- page templates generally extend `base.html`
- navigation includes Home, About, Protein Query, PROTACability, Ligand Indexer, Ligand Comparison, optional Drug GPT, and an external PROTAC Builder link
- footer copy also reinforces the companion-tool ecosystem

## Frontend and Generated Static Areas

Observed generated or runtime-managed static areas include:

- `static/charts/`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- `static/ligand_images/`

Keep these treated as generated artifacts rather than source-controlled assets.

## Codebase Conventions Observed

- `app.py` remains route-heavy and fairly monolithic
- table names in code are historically inconsistent in capitalization, but SQLite currently tolerates this
- several workflows depend on runtime-generated directories
- legacy PROTAC Builder endpoints now act as compatibility surfaces for an external companion tool rather than an in-repo builder module

## Documentation Guidance For Developers

- update public docs when route names, feature flags, or backend assumptions change
- keep PROTACability wording explicitly heuristic in UI copy and docs
- distinguish in-repo features from companion-tool handoff
- keep local databases, caches, and model artifacts out of commits

## Suggested Future Test Coverage

- route smoke tests for major public pages
- JSON contract tests for ligand lookup and PROTACability APIs
- export workflow tests using a tiny fixture database
- template rendering tests for enabled and disabled Drug GPT states
- regression tests around coordinate-cache and ligand-instance helper routes
