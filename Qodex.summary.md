# Qodex.summary

## Task
Fix PROTACability filters/results in strict RANDY backend mode.

## Original Goal
Make the PROTACability page usable on Heroku, including filters and score rendering, without pulling huge database payloads through the Heroku dyno.

## Assumptions
- Heroku V-LiSEMOD runs with `VLISMOD_DATA_BACKEND=randy`.
- Heroku points to `VLISMOD_BACKUP_URL=https://randy.rove-vernier.ts.net/backup/vlismod`.
- The RANDY token is supplied only via environment variables and must never be committed or printed.
- The production database lives on RANDY at `/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`.
- Heroku may still generate local artifacts such as CSV buffers, SVGs, and `.pml` files, but it must not depend on a production SQLite database file.
- Normal PROTACability UI requests should stay comfortably below Heroku’s request and memory limits.

## Files Inspected
- `app.py`
- `RANDY/vlismod_data_routes.py`
- `templates/protacability_assessment.html`
- `TOOLS/check_vlismod_randy_routes.py`
- `README.md`
- `docs/DEPLOYMENT.md`
- Remote runtime file: `/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`

## Files Changed
- `app.py`
- `RANDY/vlismod_data_routes.py`
- `README.md`
- `docs/DEPLOYMENT.md`
- `TOOLS/check_vlismod_randy_routes.py`
- `Qodex.summary.md`

## Files Created
- None

## Implementation Summary
- Added a response-size guard around RANDY JSON calls in `app.py` so normal UI routes can reject unexpectedly huge payloads before they become the new default path.
- Added compact RANDY PROTACability endpoints in `RANDY/vlismod_data_routes.py`:
  - `GET /backup/vlismod/protacability/filter-options`
  - `GET /backup/vlismod/protacability/search`
  - `GET /backup/vlismod/protacability/detail/<pdb_code>/<chain_id>`
  - `GET /backup/vlismod/protacability/structure-detail/<pdb_code>`
  - `GET /backup/vlismod/protacability/protein-detail`
  - `GET /backup/vlismod/protacability/target-detail`
  - `GET /backup/vlismod/protacability/export-filtered`
- Updated the Heroku-facing PROTACability routes in `app.py` so strict `randy` mode proxies those compact endpoints instead of pulling the full `/protacability/source` payload.
- Added `page/page_size` compatibility to the PROTACability search helpers while preserving the app’s existing `limit/offset` behavior.
- Synced the updated `RANDY/vlismod_data_routes.py` into the live RANDY backup receiver and reloaded Gunicorn workers.

## Key Decisions
- The full `/protacability/source` payload is not suitable for page-load filter metadata or interactive search because it can exceed 250 MB and trigger both Heroku H12 timeouts and R14/R15 memory failures.
- Filtering, sorting, grouping, pagination, and targeted detail assembly now belong on RANDY so only UI-sized JSON crosses the network to Heroku.
- Strict `randy` mode does not silently fall back to local SQLite for PROTACability routes.
- Explicit export can still be large; that is acceptable because export is a deliberate user action rather than a page-load dependency.
- Token handling stayed environment-only throughout validation.

## Commands Run
- `python -m py_compile app.py RANDY/app.py RANDY/vlismod_data_routes.py`
- `grep` and `sed` audits over PROTACability routes in `app.py`, `RANDY/vlismod_data_routes.py`, and `templates/protacability_assessment.html`
- Local Flask-client validation for RANDY blueprint responses
- Local Flask-client validation for strict `VLISMOD_DATA_BACKEND=randy` app routes
- `scp RANDY/vlismod_data_routes.py randy:/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`
- `ssh randy 'kill -HUP <gunicorn master pid>'`
- Direct `curl` validation against `https://randy.rove-vernier.ts.net/backup/vlismod`

## Validation Results
- Syntax:
  - `python -m py_compile app.py RANDY/app.py RANDY/vlismod_data_routes.py` passed.

- Direct RANDY filter-options:
  - `GET /backup/vlismod/protacability/filter-options` -> `200`
  - size: about `20,546` bytes
  - time: about `5.37s`

- Direct RANDY search:
  - `GET /backup/vlismod/protacability/search?page=1&page_size=25` -> `200`
  - size: about `88,565` bytes
  - time: about `5.42s`
  - returned `25` rows with `limit=25`, `offset=0`, and `has_more=true`

- Local strict `randy` app validation against live RANDY:
  - `GET /api/protacability/filter_options` -> `200`
  - `GET /api/protacability/search?page=1&page_size=25` -> `200`
  - `GET /api/protacability/detail/6VX2/A` -> `200`
  - `GET /api/protacability/structure_detail/6VX2?...` -> `200`
  - `GET /api/protacability/protein_detail?...` -> `200`
  - `GET /api/protacability/target_detail?...` -> `200`

- Detail payload sizes observed through the strict `randy` proxy stayed well below the old bulk source response:
  - chain detail about `8 KB`
  - structure detail about `12 KB`
  - protein detail about `39 KB`
  - target detail about `932 KB`

- Explicit filtered export remains intentionally large:
  - `GET /backup/vlismod/protacability/export-filtered?...` returned about `168 MB` for an unbounded filtered export request
  - this is acceptable for an explicit export path but not for normal UI traffic

## Known Issues
- `protacability/export-filtered` can still be very large when the filter set is broad. It is no longer used for page-load UI behavior, but it may need separate export-specific optimization later.
- The live Heroku app itself was not redeployed from this workspace during this turn, so the confirmed production fix is at the code level plus live RANDY service level. A post-deploy Heroku smoke test is still required.
- Normal PROTACability UI traffic is fixed, but export remains the one intentionally heavyweight route in this group.

## Manual Verification
1. Deploy the updated app code to Heroku.
2. Verify compact RANDY endpoints directly:
   - `curl -sS -w "\nHTTP %{http_code} bytes %{size_download} time %{time_total}\n" -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/protacability/filter-options" -o /tmp/protacability_filter_options.json`
   - `curl -sS -w "\nHTTP %{http_code} bytes %{size_download} time %{time_total}\n" -H "Authorization: Bearer $RANDY_API_TOKEN" "$VLISMOD_BACKUP_URL/protacability/search?page=1&page_size=25" -o /tmp/protacability_search.json`
3. Verify Heroku endpoints:
   - `curl -sS -w "\nHTTP %{http_code} bytes %{size_download} time %{time_total}\n" "https://vlisemod-0e358c20a94d.herokuapp.com/api/protacability/filter_options" -o /tmp/heroku_filter_options.json`
   - `curl -sS -w "\nHTTP %{http_code} bytes %{size_download} time %{time_total}\n" "https://vlisemod-0e358c20a94d.herokuapp.com/api/protacability/search?page=1&page_size=25" -o /tmp/heroku_protac_search.json`
4. Browser smoke test on `/protacability_page`:
   - confirm filter dropdowns populate
   - confirm initial results load
   - change virus/protein filters and verify results update
   - change page size and verify pagination updates
   - open structure, chain, protein, and target detail panels
   - confirm export runs only when clicked
5. Watch Heroku logs and confirm there are no new:
   - `H12 Request timeout`
   - `sqlite3.OperationalError`
   - `Error R14`
   - `Error R15`

## Suggested Next Prompt
If export performance still needs work after deploy, use:

“Optimize the PROTACability export path so large filtered exports stream or page from RANDY instead of materializing huge JSON responses, while keeping normal UI routes unchanged.”
