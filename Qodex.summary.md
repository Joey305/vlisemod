# Qodex.summary

## Task
Validate V-LiSEMOD against the real RANDY backup endpoint pattern.

## Original Goal
The user wants V-LiSEMOD to use RANDY the same way existing tools use:
- `https://randy.rove-vernier.ts.net/backup/e3`
- `https://randy.rove-vernier.ts.net/backup/hunter-job`
- `https://randy.rove-vernier.ts.net/backup/protac-event`

The goal is to confirm or implement:
- `https://randy.rove-vernier.ts.net/backup/vlismod`

before GitHub/Heroku deployment.

## Assumptions
- RANDY and V-LiSEMOD are separate Flask apps that may be deployed independently.
- The live RANDY host is reachable from this environment over HTTPS/Tailscale-style routing.
- The intended production database path on RANDY is `/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`.
- The current live RANDY service may not yet be running the latest V-LiSEMOD route prefix or env configuration.
- `VLISMOD_BACKUP_URL` should be treated as the preferred production env var, with `RANDY_API_BASE_URL` kept for legacy/local compatibility.

## Files Inspected
- `app.py`
- `RANDY/app.py`
- `RANDY/e3_data_routes.py`
- `RANDY/vlismod_data_routes.py`
- `README.md`
- `docs/DEPLOYMENT.md`
- `requirements.txt`
- `Procfile`
- `Procfile.randy`
- `.gitignore`

## Files Changed
- `app.py`
- `RANDY/app.py`
- `RANDY/vlismod_data_routes.py`
- `README.md`
- `docs/DEPLOYMENT.md`

## Files Created
- `Qodex.summary.md`

## Implementation Summary
Inspected the existing RANDY route conventions and confirmed that E3 already uses a backup-style blueprint under `/backup/e3`, while V-LiSEMOD had only been exposed under `/api/vlismod`.

Updated the RANDY V-LiSEMOD route module so it now exposes both:
- `/api/vlismod/*`
- `/backup/vlismod/*`

Updated the main V-LiSEMOD client helper so production can use:
- `VLISMOD_BACKUP_URL=https://randy.rove-vernier.ts.net/backup/vlismod`

while preserving:
- `RANDY_API_BASE_URL`
- `RANDY_API_TOKEN`
- `VLISMOD_DATA_BACKEND`

The client now avoids double-prefix bugs by appending route names like `viruses` directly under the configured base URL.

## Key Decisions
- Live RANDY URL testing is required because local Flask-client tests only prove code paths, not deployed routing, token configuration, or database wiring.
- Local Flask-client checks were still used as secondary validation after code changes, but not treated as final proof.
- `/api/vlismod` was preserved as a compatibility alias because the live RANDY host is already serving that prefix today.
- `VLISMOD_BACKUP_URL` should be the preferred production env var going forward.

## Commands Run
- Repo inspection commands:
  - `pwd`
  - `git status --short`
  - `ls`
  - `ls RANDY`
  - `grep` searches for backup route patterns, blueprints, env vars, and `viral_data.db`
  - `find` for `Procfile*`, `.gitignore`, and runtime files
- Code inspection commands with `sed`
- Syntax check:
  - `python -m py_compile app.py RANDY/app.py RANDY/e3_data_routes.py RANDY/vlismod_data_routes.py`
- Live HTTPS checks with `curl` and Python `requests` against:
  - `https://randy.rove-vernier.ts.net/backup/e3/healthz`
  - `https://randy.rove-vernier.ts.net/backup/vlismod/health`
  - `https://randy.rove-vernier.ts.net/api/vlismod/health`
  - `https://randy.rove-vernier.ts.net/api/vlismod/viruses`
  - `https://randy.rove-vernier.ts.net/api/vlismod/db-health`
  - `https://randy.rove-vernier.ts.net/api/vlismod/ligands/list`
- Secondary local verification with Flask `test_client()` to confirm the new `/backup/vlismod` alias exists in code.

## Validation Results
- Live RANDY URL tested:
  - `https://randy.rove-vernier.ts.net/backup/vlismod`
  - `https://randy.rove-vernier.ts.net/api/vlismod`

- Live HTTPS status results:
  - `GET /backup/e3/healthz` -> `401`
  - `GET /backup/vlismod/health` -> `404`
  - `GET /backup/vlismod/viruses` -> `404`
  - `GET /api/vlismod/health` -> `200`
  - `GET /api/vlismod/viruses` -> `500`
  - `GET /api/vlismod/db-health` -> `500`
  - `GET /api/vlismod/ligands/list` -> `500`

- Live response interpretation:
  - The live host is reachable.
  - `/backup/vlismod` does not currently exist on the deployed RANDY service.
  - `/api/vlismod` does exist on the deployed RANDY service.
  - The deployed RANDY V-LiSEMOD service is currently misconfigured:
    - `auth_configured` is `false`
    - `db_exists` is `false`
    - `db_path` is `/home/jxs794/PROTAC_BUILDER/viral_data.db`
  - Because the live service has no configured V-LiSEMOD token, live unauthorized behavior for V-LiSEMOD could not be validated as `401`; the deployed `/api/vlismod/*` routes return `500` instead.

- V-LiSEMOD strict `randy` mode using live `VLISMOD_BACKUP_URL`:
  - `/get_viruses` -> `502`
  - `/get_ligands_list` -> `502`
  - `/get_pdb_codes/<known virus>` -> `502`
  - `/get_ligands/<known pdb>` -> `502`
  - `/check_functional_groups/<known pdb>` -> `502`
  - `/get_pdb_mapping/<known ligand>` -> `502`
  - Error is now explicit: `RANDY API request failed with status 404.`

- Secondary local code validation:
  - `/api/vlismod/health` exists
  - `/backup/vlismod/health` exists
  - `/backup/vlismod/viruses` works with bearer auth when the local app is configured

- `VLISMOD_DATA_BACKEND=auto` fallback still works locally when RANDY is unavailable.

## Known Issues
- The deployed live RANDY service does not currently expose `/backup/vlismod`.
- The deployed live RANDY V-LiSEMOD service is using the wrong database path and does not have a configured V-LiSEMOD token.
- Live unauthorized `401` behavior for V-LiSEMOD cannot pass until the remote service is redeployed or reconfigured.
- Strict V-LiSEMOD `randy` mode cannot pass against the live backup URL until `/backup/vlismod` is actually present on RANDY.
- Production deployment still needs a durable strategy for the SQLite database outside ephemeral Heroku dyno storage.

## Manual Verification
If you want to verify the remote service from a machine that can reach RANDY, run:

```bash
export VLISMOD_BACKUP_URL="https://randy.rove-vernier.ts.net/backup/vlismod"
export RANDY_API_TOKEN="<token from env only>"

curl -i "$VLISMOD_BACKUP_URL/health"
curl -i "$VLISMOD_BACKUP_URL/viruses"
curl -i -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/db-health"
curl -i -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/viruses"
curl -i -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/ligands/list"
```

If the live service still only exposes `/api/vlismod`, first confirm current state:

```bash
curl -i "https://randy.rove-vernier.ts.net/api/vlismod/health"
curl -i -H "Authorization: Bearer $RANDY_API_TOKEN" "https://randy.rove-vernier.ts.net/api/vlismod/db-health"
```

On RANDY itself, also verify the runtime env:

```bash
echo "$VLISMOD_DB_PATH"
ls -lh /home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db
test -r /home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db && echo readable
```

## Suggested Next Prompt
Set up the final RANDY/V-LiSEMOD deployment alignment by updating the live RANDY service to load the new `/backup/vlismod` routes, point `VLISMOD_DB_PATH` at `/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`, configure the server-side token env var, and then re-run live HTTPS validation before preparing GitHub and Heroku deployment steps.
