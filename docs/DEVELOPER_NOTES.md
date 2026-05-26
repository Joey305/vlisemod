# Developer Notes

## Flask App Structure

The main application lives in `app.py` and exposes most public routes directly from a single Flask app instance. An optional Drug GPT assistant module lives in `DRUGapp.py` and is registered as the `dp` blueprint at `/drugapp` only when `ENABLE_DRUG_GPT` is enabled.

The public app is no-login by default. Authentication routes, `users.db`, and Flask-Login have been removed from the runtime path.

## Route Inventory Instructions

To inspect routes in a fully provisioned environment:

```bash
python -m flask --app app routes
```

or

```bash
flask --app app routes
```

During this documentation pass, route inspection via Flask CLI failed in the current environment because `seaborn` was missing, so route references were verified from source using code search.

## High-Level Route Areas

Key public routes currently include:

- `/`
- `/about`
- `/ligand_indexer`
- `/compare_ligands`
- `/query_protein_virus_page`
- `/protacability_page`
- `/healthz`

Key JSON/API areas currently include:

- ligand selection and mapping endpoints
- export endpoints
- PROTACability filter/search/detail/export endpoints
- coordinate and ligand-instance debug/viewer endpoints

## Template / Navigation Conventions

- shared shell: `templates/base.html`
- page templates extend `base.html`
- nav currently includes Home, About, Protein Query, PROTACability, Ligand Indexer, Ligand Comparison, optional Drug GPT, and an external PROTAC Builder link
- the footer also reinforces the companion-tool ecosystem

## Optional Module Handling

Current flags in `app.py`:

- `SHOW_DRUG_GPT_NAV`
- `ENABLE_DRUG_GPT`
- `ENABLE_LOCAL_LLM`

If Drug GPT is disabled, `app.py` still serves disabled-state routes for `/drugapp/` and `/drugapp/query` so the UI can fail gracefully instead of hard-404ing.

## Drug GPT / BioGPT Disabled-State Route Issue

Several templates still probe `/drugapp/` with a `HEAD` request and then show either an embedded assistant or an unavailable-state popup. That means:

- disabled does not mean absent,
- UI behavior depends on both feature flags and route availability,
- public docs should describe the module as optional and deployment-dependent.

## Light / Dark Mode Notes

`templates/base.html` sets a `data-theme` attribute from `localStorage` key `vlisemod_theme`, and shared styling/scripts support a light/dark toggle. When editing templates or shared CSS, preserve theme-aware variables and contrast-sensitive surfaces.

## Frontend / Static Organization

Observed frontend organization includes:

- `static/css/styles.css`
- `static/js/theme-toggle.js`
- `static/js/scripts.js`
- `static/js/ngl_viewer_helpers.js`

Additional static/generated areas are used for charts, coordinate caches, ligand SDF caches, and ligand images.

## Codebase Conventions Discovered During Inspection

- the codebase is route-heavy and still fairly monolithic in `app.py`
- table names in code are historically inconsistent in capitalization, but SQLite currently tolerates this
- several workflows rely on runtime-generated folders
- legacy PROTAC Builder endpoints now return deprecation payloads and redirect users to the standalone tool
- documentation should distinguish clearly between maintained in-repo modules and external companion tools

## Auth Note

`users.db` is no longer part of the application runtime. Keep any lingering documentation references limited to deprecation or migration notes only.

## Suggested Future Test Coverage

- route smoke tests for the main public pages
- JSON contract tests for ligand lookup and PROTACability APIs
- export workflow tests with a tiny fixture database
- template rendering tests for enabled and disabled Drug GPT states
- regression tests for coordinate-cache and ligand-instance lookup helpers

## Developer Workflow Advice

- use `rg` or route search before changing copy that names pages or endpoints
- update docs when route names or env flags change
- keep generated data and local-only assets out of commits
- treat PROTACability language carefully in code comments, docs, and UI copy
