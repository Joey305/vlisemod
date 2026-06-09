# Deployment Guide

## Local Development Overview

V-LiSEMOD is a public Flask application backed by `viral_data.db` and a set of locally available structure-derived assets. A minimal deployment requires the Python dependencies, the expected scientific data files, and any environment variables needed for optional features.

## Environment Variables

Confirmed from the current codebase:

- `FLASK_SECRET_KEY`: session secret; should be overridden outside local testing
- `VLISMOD_DATA_BACKEND`: data source mode for the simple lookup routes; `local`, `randy`, or `auto`
- `VLISMOD_BACKUP_URL`: preferred production-style RANDY V-LiSEMOD base URL, for example `https://randy.rove-vernier.ts.net/backup/vlismod`
- `RANDY_API_BASE_URL`: base URL for the RANDY service, for example `http://127.0.0.1:5001`
- `RANDY_API_TOKEN`: bearer token used by V-LiSEMOD when calling RANDY
- `SHOW_DRUG_GPT_NAV`: show or hide the Drug GPT navigation entry
- `ENABLE_DRUG_GPT`: enable or disable the `/drugapp/` blueprint
- `ENABLE_LOCAL_LLM`: allow local model loading for the optional assistant module
- `MODEL_ID` or `LLM_MODEL_ID`: local assistant model identifier
- `PROTAC_BUILDER_EXTERNAL_URL`: override default external handoff URL

Optional Drug GPT-related variables in `DRUGapp.py` also include:

- `HF_TOKEN`, `HUGGINGFACE_HUB_TOKEN`, or `HUGGING_FACE_HUB_TOKEN`
- `GEN_MAX_CONCURRENT`
- `GEN_ACQUIRE_TIMEOUT`

Additional variables used by the RANDY backend:

- `VLISMOD_API_TOKEN`: bearer token expected by RANDY for `/api/vlismod/*`
- `VLISMOD_DB_PATH`: filesystem path to the SQLite database RANDY should query
- `RANDY_BACKUP_TOKEN` or `PROTAC_BACKUP_TOKEN`: accepted as fallback token env vars for `/backup/vlismod/*` if you want V-LiSEMOD to follow the same backup-token convention as the other RANDY services

## Running With `python app.py`

If the codebase exposes the usual Flask `app` object and the runtime dependencies are installed, a basic local start path is:

```bash
python app.py
```

Because the repository uses direct imports such as `matplotlib`, `seaborn`, `rdkit`, `pandas`, and optional `biopython`, missing packages will prevent startup.

Current source note:

- the checked-in `if __name__ == "__main__"` block runs the dev server on `127.0.0.1:5003`

## Running With Waitress or Gunicorn

The top of `app.py` includes a note showing a Waitress-style entrypoint:

```bash
waitress-serve --listen=127.0.0.1:5002 app:app
```

This means the documented Waitress example and the current `python app.py` dev port are not identical; choose one intentionally for your environment.

No `Procfile` is currently present in the repository root, and `gunicorn` is not listed in `requirements.txt`. If you want Gunicorn in production, add and test it explicitly rather than assuming it is already configured.

## Procfile Notes

- The repository root now includes `Procfile` for the main V-LiSEMOD Flask app using `gunicorn app:app`.
- A separate `Procfile.randy` is also available for RANDY-only deployment and uses `gunicorn app:APP --chdir RANDY ...`.
- Heroku only consumes a root `Procfile`, so separate Heroku apps still need an explicit deployment workflow:
  - deploy V-LiSEMOD with the root `Procfile`
  - deploy RANDY from a repo layout, branch, or build step that makes the `Procfile.randy` command the active web process

## GitHub-First Deployment Workflow

Recommended public-repo workflow:

1. Keep source, templates, and docs in GitHub.
2. Keep `viral_data.db`, generated caches, model weights, and large structure assets out of the repository.
3. Provision the runtime environment first.
4. Upload or mount required data artifacts separately.
5. Set environment variables in the hosting platform.
6. Validate `/healthz` and major public pages after deploy.

## Heroku Dashboard Deployment From GitHub

There is no Heroku-specific config file in the inspected root, but the app can still be deployed from GitHub through the Heroku dashboard if:

1. the Python buildpack is configured,
2. a correct start command is provided,
3. required large data assets are provisioned outside Git,
4. environment variables are configured in the dashboard, and
5. writable directories exist for generated outputs if the platform filesystem allows them.

Do not commit Heroku tokens, CLI auth artifacts, or platform-specific secrets.

### RANDY as a Separate Service

If RANDY is deployed separately from V-LiSEMOD, treat them as two Flask services with coordinated environment variables.

RANDY service:

- set `VLISMOD_API_TOKEN`
- set `VLISMOD_DB_PATH`
- confirm `GET /backup/vlismod/health` and authenticated `GET /backup/vlismod/db-health`
- preserve `/api/vlismod/*` only as a compatibility alias if useful

V-LiSEMOD service:

- set `VLISMOD_DATA_BACKEND=randy` for strict remote mode or `VLISMOD_DATA_BACKEND=auto` for fallback mode
- set `VLISMOD_BACKUP_URL=https://randy.rove-vernier.ts.net/backup/vlismod`
- set `RANDY_API_TOKEN`
- optionally keep `RANDY_API_BASE_URL` only for local/dev compatibility
- keep `viral_data.db` available only if you want local fallback when using `auto`

Suggested local verification flow:

```bash
export VLISMOD_API_TOKEN="dev-token"
export VLISMOD_DB_PATH="../viral_data.db"
python RANDY/app.py
```

Then in a second shell:

```bash
export VLISMOD_DATA_BACKEND=auto
export VLISMOD_BACKUP_URL="http://127.0.0.1:8787/backup/vlismod"
export RANDY_API_TOKEN="dev-token"
python app.py
```

Base-URL precedence in the app:

1. `VLISMOD_BACKUP_URL`
2. `RANDY_API_BASE_URL`
3. local-only behavior when RANDY mode is not configured

If you set `VLISMOD_BACKUP_URL`, the client appends route names like `viruses` or `pdb-codes` directly under that base, which avoids `/api/vlismod` double-prefix bugs.

### Storage Caveat For Heroku

RANDY currently reads a SQLite file from `VLISMOD_DB_PATH`. That is workable for local testing and for a persistent non-Heroku host, but it is not a durable Heroku storage strategy by itself.

Key implications:

- Heroku dyno filesystems are ephemeral, so a database copied onto a dyno will not be durable across restarts or deploys.
- Bundling a large SQLite file into the slug risks slug-size and startup-time issues, and the file remains effectively immutable at runtime.
- If the durable database stays on another machine, Heroku will need a reachable API layer or another transport path to access that data safely.
- If long-term hosted durability is required, plan a persistent database or object-storage-backed strategy before production rollout.

## Azure Notes

No Azure-specific deployment configuration was identified during inspection. If Azure App Service or a container-based deployment is used later, keep in mind:

- the app depends on local writable/cache directories unless adapted,
- SQLite and large structure assets may not be ideal for highly ephemeral instances,
- optional local-model serving is likely too heavy for default lightweight deployments.

## Health Check Route

Confirmed route:

- `/healthz`

Use this as a basic deployment smoke check after startup.

## Large File / Data Handling

The repository expects significant local data outside the public source footprint:

- `viral_data.db`
- `PDB_FILES/`
- generated caches and exports in `static/` and output folders

Treat these as provisioned runtime assets. Public GitHub should carry code and docs, not bulk regenerated datasets.

## Common Deployment Pitfalls

### Missing database

- Symptom: query pages load poorly or endpoints fail when `viral_data.db` is absent or incomplete.
- Fix: provide a valid local database before starting the app.

### Missing environment variables

- Symptom: weak anonymous session security, incorrect external handoff URL, or optional modules not behaving as intended.
- Fix: set `FLASK_SECRET_KEY` and any feature flags explicitly in the deployment environment.

### Optional model loading

- Symptom: `/drugapp/` fails or local startup becomes heavy.
- Fix: keep `ENABLE_DRUG_GPT=0` and `ENABLE_LOCAL_LLM=0` for lightweight public deployments unless the model runtime is intentionally provisioned.

### Drug GPT disabled/enabled routes

- Symptom: UI surfaces check `/drugapp/` and may show an unavailable assistant state.
- Fix: this is expected when the module is disabled. If enabling it, verify both the nav flags and runtime dependencies.

### Generated folders not writable

- Symptom: chart generation, cached coordinates, ligand SDF generation, or session exports fail.
- Fix: ensure runtime write access for generated output directories or refactor paths for your hosting platform.
