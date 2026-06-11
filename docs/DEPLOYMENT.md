# Deployment Guide

## Overview

V-LiSEMOD is a Flask application that depends on provisioned scientific data, writable output locations, and optional environment flags for deployment-specific behavior. Public source availability does not mean the application is self-contained without data provisioning.

## Confirmed Environment Variables

From the current codebase:

- `FLASK_SECRET_KEY`
- `VLISMOD_LOCAL_DB_PATH`
- `VLISMOD_DATA_BACKEND`
- `VLISMOD_BACKUP_URL`
- `RANDY_API_BASE_URL`
- `RANDY_API_TOKEN`
- `RANDY_API_TIMEOUT_SECONDS`
- `SHOW_DRUG_GPT_NAV`
- `ENABLE_DRUG_GPT`
- `ENABLE_LOCAL_LLM`
- `MODEL_ID`
- `LLM_MODEL_ID`
- `PROTAC_BUILDER_EXTERNAL_URL`

Optional Drug GPT-related variables imported in `DRUGapp.py`:

- `HF_TOKEN`
- `HUGGINGFACE_HUB_TOKEN`
- `HUGGING_FACE_HUB_TOKEN`
- `GEN_MAX_CONCURRENT`
- `GEN_ACQUIRE_TIMEOUT`

Additional variables used by the RANDY service:

- `VLISMOD_API_TOKEN`
- `VLISMOD_DB_PATH`
- `RANDY_BACKUP_TOKEN`
- `PROTAC_BACKUP_TOKEN`

## Local Development Run

A minimal local run path is:

```bash
python app.py
```

The current source uses the Flask app object in `app.py` and the checked-in dev server configuration binds to `127.0.0.1:5003` when run directly.

## Alternative Entrypoints

The source also contains a Waitress-style example near the top of `app.py`:

```bash
waitress-serve --listen=127.0.0.1:5002 app:app
```

The repository root includes a `Procfile` using `gunicorn app:app --timeout 900`. Treat production entrypoints as deployment choices that still need environment-specific validation.

## Local SQLite Mode

Local mode uses `viral_data.db` directly. At minimum, confirm:

1. the database file exists,
2. the expected tables are present,
3. generated output folders are writable, and
4. optional modules remain disabled unless intentionally provisioned.

## RANDY-Backed Mode

Selected route groups can run against a separate RANDY service using `VLISMOD_DATA_BACKEND=randy` or `auto`.

Operationally:

- `VLISMOD_BACKUP_URL` is preferred when available,
- `RANDY_API_BASE_URL` remains useful for local development compatibility,
- `RANDY_API_TOKEN` or compatible fallback tokens must be supplied,
- normal PROTACability UI flows should use the compact filter, search, and detail endpoints rather than a bulk source payload.

Keep the documentation phrased as deployment options, not as a claim that every public deployment includes the same backend topology.

## Health and Smoke Checks

Confirmed route:

- `/healthz`

Recommended smoke checks after startup:

1. `/healthz`
2. `/`
3. `/query_protein_virus_page`
4. `/ligand_indexer`
5. `/compare_ligands`
6. `/protacability_page`
7. `/drugapp/` only if the optional module is intentionally enabled

## Data Provisioning Expectations

The application expects significant local or provisioned data outside the public source tree, including:

- `viral_data.db`
- `PDB_FILES/`
- generated caches and export directories

Treat these as runtime assets. Public GitHub should carry source and documentation, not local databases or bulk regenerated outputs.

## Common Deployment Pitfalls

### Missing database

- Symptom: query pages or API routes return empty or degraded results.
- Fix: provide a valid local database or a working RANDY-backed data source.

### Missing generated-folder permissions

- Symptom: chart generation, ligand imagery, cached coordinates, or PyMOL-oriented exports fail.
- Fix: ensure writable directories exist for runtime-generated artifacts.

### Optional assistant accidentally enabled

- Symptom: startup becomes heavy or `/drugapp/` fails due to missing model/runtime dependencies.
- Fix: keep `ENABLE_DRUG_GPT=0` and `ENABLE_LOCAL_LLM=0` for lightweight public deployments unless the model runtime is intentionally provisioned.

### RANDY mode misconfiguration

- Symptom: lookup and PROTACability routes fail or time out.
- Fix: verify base URL, token configuration, timeout, and route compatibility between V-LiSEMOD and RANDY.

## Production-Readiness Caveat

This repository supports research-oriented deployment, but public documentation should not imply production-hardened infrastructure, durable hosted storage, formal API guarantees, or enterprise security controls unless those have been implemented and separately validated.
