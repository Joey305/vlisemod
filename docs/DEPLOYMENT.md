# Deployment Guide

## Local Development Overview

V-LiSEMOD is a public Flask application backed by `viral_data.db` and a set of locally available structure-derived assets. A minimal deployment requires the Python dependencies, the expected scientific data files, and any environment variables needed for optional features.

## Environment Variables

Confirmed from the current codebase:

- `FLASK_SECRET_KEY`: session secret; should be overridden outside local testing
- `SHOW_DRUG_GPT_NAV`: show or hide the Drug GPT navigation entry
- `ENABLE_DRUG_GPT`: enable or disable the `/drugapp/` blueprint
- `ENABLE_LOCAL_LLM`: allow local model loading for the optional assistant module
- `MODEL_ID` or `LLM_MODEL_ID`: local assistant model identifier
- `PROTAC_BUILDER_EXTERNAL_URL`: override default external handoff URL

Optional Drug GPT-related variables in `DRUGapp.py` also include:

- `HF_TOKEN`, `HUGGINGFACE_HUB_TOKEN`, or `HUGGING_FACE_HUB_TOKEN`
- `GEN_MAX_CONCURRENT`
- `GEN_ACQUIRE_TIMEOUT`

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

- No root `Procfile` was found during inspection.
- If deploying to a platform that expects a `Procfile`, add one intentionally and keep it aligned with the true runtime command.

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
