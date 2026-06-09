# Qodex.summary

## Task
Test live RANDY `protac-backup` V-LiSEMOD endpoints after restart.

## Original Goal
Confirm that after copying the files to RANDY and restarting `protac-backup`, V-LiSEMOD works through `https://randy.rove-vernier.ts.net/backup/vlismod`, then decide whether it is safe to proceed to GitHub/Heroku setup.

## Assumptions
- `protac-backup.service` is the live Gunicorn service for RANDY and serves `backup_receiver.app:APP`.
- The intended production database path is `/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`.
- `VLISMOD_API_TOKEN` and `VLISMOD_DB_PATH` should be loaded from `/home/jxs794/PROTAC_BUILDER/.env`.
- `/backup/vlismod` is the intended production route prefix.
- The V-LiSEMOD app source itself is not present on RANDY; only the backup receiver and the database directory are present there.

## Files Inspected
- Remote: `/home/jxs794/PROTAC_BUILDER/backup_receiver/app.py`
- Remote: `/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`
- Remote: `/etc/systemd/system/protac-backup.service`
- Remote: `/home/jxs794/PROTAC_BUILDER/.env`
- Remote: `/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`
- Local: `app.py`
- Local: `.gitignore`
- Local: `requirements.txt`
- Local: `Procfile`
- Local: `Procfile.randy`

## Files Changed
- Remote runtime config: `/home/jxs794/PROTAC_BUILDER/.env`
- Local summary: `Qodex.summary.md`

## Files Created
- None

## Implementation Summary
This was mostly a live verification pass with one runtime configuration fix on RANDY. The deployed `backup_receiver` code already included the new `/backup/vlismod` blueprint and the `/api/vlismod` compatibility alias. The live failure after restart was caused by the service still reading the wrong default database path.

I updated the remote `.env` to set:
- `VLISMOD_DB_PATH=/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`
- `VLISMOD_API_TOKEN=<set>`

Because `sudo` restart was not available non-interactively, I sent `HUP` to the Gunicorn master process owned by `jxs794`, which successfully booted fresh workers and picked up the `.env` changes.

## Key Decisions
- `/backup/vlismod` is now the correct production URL.
- Live HTTPS validation was treated as the source of truth; local Flask-client checks were only secondary.
- Heroku/GitHub setup can proceed from an integration standpoint because the live RANDY endpoint now works, but deployment still depends on the usual secret-handling and persistent-database considerations.

## Commands Run
- Remote service inspection:
  - `ssh randy 'pwd; hostname; ...'`
  - `ssh randy 'systemctl status protac-backup --no-pager'`
  - `ssh randy 'systemctl cat protac-backup'`
  - `ssh randy 'systemctl show protac-backup -p Environment --no-pager'`
  - `ssh randy 'ps -ef | grep gunicorn | grep -v grep'`
- Remote code inspection with `sed` and `grep`
- Remote DB checks with `/home/jxs794/PROTAC_BUILDER/.venv/bin/python` using SQLite
- Remote internal endpoint checks with `curl` against `http://127.0.0.1:8787/backup/vlismod`
- Remote external endpoint checks with `curl` against `https://randy.rove-vernier.ts.net/backup/vlismod`
- Remote `.env` update with a short Python script
- Gunicorn worker reload with `kill -HUP <master-pid>`
- Local strict `randy`-mode V-LiSEMOD route checks against the live RANDY URL
- Local git safety checks

## Validation Results
- `protac-backup.service`: active/running
- Gunicorn target: `backup_receiver.app:APP`
- Gunicorn bind: `127.0.0.1:8787`

- Internal localhost checks:
  - `/backup/vlismod/health` -> `200`
  - `/backup/vlismod/viruses` without token -> `401`
  - `/backup/vlismod/db-health` with token -> `200`
  - `/backup/vlismod/viruses` with token -> `200`
  - `/backup/vlismod/ligands/list` with token -> `200`

- External HTTPS checks:
  - `/backup/vlismod/health` -> `200`
  - `/backup/vlismod/viruses` without token -> `401`
  - `/backup/vlismod/db-health` with token -> `200`
  - `/backup/vlismod/viruses` with token -> `200`
  - `/backup/vlismod/ligands/list` with token -> `200`
  - `/backup/vlismod/pdb-codes?virus_name=Human%20immunodeficiency%20virus%201` -> `200`
  - `/backup/vlismod/ligands?pdb_code=1A8G` -> `200`
  - `/backup/vlismod/functional-groups/check?pdb_code=1R6N` -> `200`
  - `/backup/vlismod/pdb-mapping?ligand_code=2Z4` -> `200`

- `db-health` now confirms:
  - `db_exists=true`
  - `db_path=/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`
  - `auth_configured=true` via `/health`

- V-LiSEMOD strict `randy` mode against live RANDY:
  - `/get_viruses` -> `200`
  - `/get_ligands_list` -> `200`
  - `/get_pdb_codes/Human immunodeficiency virus 1` -> `200`
  - `/get_ligands/1A8G` -> `200`
  - `/check_functional_groups/1R6N` -> `200`
  - `/get_pdb_mapping/2Z4` -> `200`

- Local fallback:
  - Not used during strict `randy` mode verification.

## Known Issues
- `systemctl restart protac-backup` could not be run non-interactively through `sudo`; worker reload was done with `HUP` instead.
- The actual V-LiSEMOD app source tree is not present on RANDY, so the strict `randy`-mode app check was performed from the local V-LiSEMOD codebase against the live RANDY URL rather than from a remote app checkout.
- The repository still has broader pre-existing local changes unrelated to this validation flow; review before pushing.

## Manual Verification
From RANDY:

```bash
cd /home/jxs794/PROTAC_BUILDER
grep -E "^(VLISMOD_API_TOKEN|VLISMOD_DB_PATH)=" .env | sed 's/=.*/=<set>/'

TOKEN="<token from env only>"
BASE="http://127.0.0.1:8787/backup/vlismod"

curl -i "$BASE/health"
curl -i "$BASE/viruses"
curl -i -H "Authorization: Bearer $TOKEN" "$BASE/db-health"
curl -i -H "Authorization: Bearer $TOKEN" "$BASE/viruses"
curl -i -H "Authorization: Bearer $TOKEN" "$BASE/ligands/list"
```

From any machine that can reach the Tailscale URL:

```bash
TOKEN="<token from env only>"
BASE="https://randy.rove-vernier.ts.net/backup/vlismod"

curl -i "$BASE/health"
curl -i "$BASE/viruses"
curl -i -H "Authorization: Bearer $TOKEN" "$BASE/db-health"
curl -i -H "Authorization: Bearer $TOKEN" "$BASE/viruses"
curl -i -H "Authorization: Bearer $TOKEN" "$BASE/ligands/list"
```

## Suggested Next Prompt
Prepare the GitHub and Heroku deployment steps for V-LiSEMOD now that the live RANDY backup endpoint is working, including a clean review of staged files, final env-var documentation, and the exact Heroku config vars needed for `VLISMOD_DATA_BACKEND=randy` with `VLISMOD_BACKUP_URL=https://randy.rove-vernier.ts.net/backup/vlismod`.
