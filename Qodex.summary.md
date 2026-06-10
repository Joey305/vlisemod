# Qodex.summary

## Task
Fix PROTACability metric clarity, tab behavior, and loading UX.

## Original Goal
The user wants PROTACability filters/results to feel responsive and trustworthy, wants to verify whether “Mapped/exposed warheads = 6” is correct, and wants a themed loading wheel with rotating messages during slow page/data loads.

## Assumptions
- Heroku V-LiSEMOD runs with `VLISMOD_DATA_BACKEND=randy`.
- Heroku points to `VLISMOD_BACKUP_URL=https://randy.rove-vernier.ts.net/backup/vlismod`.
- The RANDY token is supplied only through environment variables and was never printed or committed.
- The production database lives on RANDY at `/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`.
- PROTACability grouped views (`targets`, `protein`, `summary`) are intended to summarize current filtered result groups rather than raw warhead-table totals.
- The PROTACability page should stay compact-RANDY-backed and must not reintroduce the huge `/protacability/source` payload into normal UI flows.

## Files Inspected
- `/Users/jxs794/Documents/VLISEMOD/app.py`
- `/Users/jxs794/Documents/VLISEMOD/RANDY/vlismod_data_routes.py`
- `/Users/jxs794/Documents/VLISEMOD/templates/protacability_assessment.html`
- `/Users/jxs794/Documents/VLISEMOD/static/css/styles.css`
- `/Users/jxs794/Documents/VLISEMOD/Qodex.summary.md`
- Remote runtime file: `/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`

## Files Changed
- `/Users/jxs794/Documents/VLISEMOD/app.py`
- `/Users/jxs794/Documents/VLISEMOD/RANDY/vlismod_data_routes.py`
- `/Users/jxs794/Documents/VLISEMOD/templates/protacability_assessment.html`
- `/Users/jxs794/Documents/VLISEMOD/static/css/styles.css`
- `/Users/jxs794/Documents/VLISEMOD/Qodex.summary.md`

## Files Created
- None

## Implementation Summary
- Audited the `candidate_warheads_with_exposed_mapped_atoms` summary card and confirmed the old `6` value came from counting only grouped representative rows in `targets` and `protein` views.
- Added explicit grouped evidence flags in both `app.py` and `RANDY/vlismod_data_routes.py` so grouped `targets` and `protein` summaries now count any grouped structure context with mapped + solvent-exposed mapped atoms.
- Left `summary` and `chains` semantics intact, because structure-group and raw chain-row counts already matched the intended evidence scope.
- Reworded the second summary-card label per view so users can tell whether they are looking at target groups, protein groups, structures, or chain rows.
- Added a reusable full-page PROTACability loading overlay with rotating messages, reduced-motion handling, light/dark theme support, and error-safe hide behavior.
- Wrapped initial bootstrap, filter refresh, view switching, search, pagination, reset/apply flows, and detail fetches in loading helpers.
- Added lightweight backend diagnostics for filter/search payload size, elapsed time, totals, and mapped/exposed summary counts.
- Fixed a local fallback bug in `/api/protacability/target_detail` where the code referenced `conn` outside a DB context.
- Synced the updated RANDY route file into `/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py` and reloaded Gunicorn.

## Metric Audit Summary
- Original displayed metric:
  - `Mapped/exposed warheads = 6`
- Direct SQL comparison on RANDY `protacability_warhead_linkability`:
  - total warhead rows: `92482`
  - mapped rows: `3468`
  - mapped + solvent-exposed rows: `2659`
  - distinct structures with mapped + solvent-exposed evidence: `2513`
  - distinct target groups with mapped + solvent-exposed evidence: `61`
- Local grouped-result audit:
  - old grouped target/protein summary logic: `6`
  - grouped target/protein rows with any qualifying grouped evidence: `20`
  - grouped structure summary count: `2509`
  - raw chain-row count: `8935`
- Final metric definition:
  - `targets`: grouped target rows in the current filtered result set with at least one mapped + solvent-exposed warhead-supporting structure context
  - `protein`: grouped protein rows in the current filtered result set with at least one mapped + solvent-exposed warhead-supporting structure context
  - `summary`: grouped structure rows with mapped + solvent-exposed mapped atoms
  - `chains`: raw chain rows with mapped + solvent-exposed mapped atoms
- Final label and helper text:
  - `targets`: `Mapped/exposed target groups`
  - `protein`: `Mapped/exposed protein groups`
  - `summary`: `Mapped/exposed structures`
  - `chains`: `Mapped/exposed chain rows`

## Key Decisions
- `6` was not a trustworthy “total mapped/exposed warheads” number. It undercounted grouped target/protein evidence because it depended on whichever representative row won the grouping priority sort.
- The backend summary needed a real fix for grouped `targets` and `protein` views, not just a wording tweak.
- The final labels were still narrowed to per-view language so users do not mistake grouped evidence counts for raw database warhead totals.
- The loading overlay is page-level because initial bootstrap, grouped search, and tab switches all replace the main evidence view and can otherwise look stalled.
- Theme compatibility was handled with existing site variables (`--page-bg`, `--card-bg`, `--accent-green`, `--accent-orange`, `--text-muted`) so the loader fits both light and dark mode without a separate palette.

## Commands Run
- `pwd`
- `git status --short`
- `python -m py_compile app.py RANDY/app.py RANDY/vlismod_data_routes.py`
- `grep` audits over PROTACability metrics, search, and view handlers in `app.py`, `RANDY/vlismod_data_routes.py`, `templates/protacability_assessment.html`, and `static/css/styles.css`
- Local Python checks against `_prepare_protacability_result_set_from_rows(...)`
- Local Flask test-client checks for PROTACability routes
- `node --check /tmp/protacability_assessment.js`
- Remote RANDY SQL checks over `protacability_warhead_linkability`
- `scp RANDY/vlismod_data_routes.py randy:/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`
- `ssh randy 'kill -HUP 4020670'`
- Remote Python validation inside `/home/jxs794/PROTAC_BUILDER/backup_receiver`

## Validation Results
- Syntax:
  - `python -m py_compile app.py RANDY/app.py RANDY/vlismod_data_routes.py` passed.
- Template JS syntax:
  - extracted the inline script to `/tmp/protacability_assessment.js`
  - `node --check /tmp/protacability_assessment.js` passed.
- Metric validation:
  - grouped target/protein summary count now returns `20` instead of `6`
  - grouped structure summary remains `2509`
  - chain-row summary remains `8935`
- Local route validation:
  - `GET /api/protacability/filter_options` -> `200`
  - `GET /api/protacability/search?view=targets&page=1&page_size=50` -> `200`
  - `GET /api/protacability/search?view=protein&page=1&page_size=50` -> `200`
  - `GET /api/protacability/search?view=summary&page=1&page_size=50` -> `200`
  - `GET /api/protacability/search?view=chains&page=1&page_size=50` -> `200`
  - `GET /api/protacability/detail/6VX2/A` -> `200`
  - `GET /api/protacability/protein_detail?...` -> `200`
  - `GET /api/protacability/target_detail?...` -> `200`
- View behavior validation:
  - all four views still send and return distinct `view` values (`targets`, `protein`, `summary`, `chains`)
  - each view now returns distinct summary values and row sets
  - sort-reset and offset-reset behavior still occurs through `selectView(...)` + search
- Loading overlay behavior:
  - overlay is now wired to initial bootstrap, filter refresh, apply/reset, view switching, pagination, search, and detail fetches
  - rotating messages are configured at `4s`
  - reduced-motion CSS disables spinner animation and transitions
- Live RANDY runtime validation:
  - remote grouped summary calculation inside `/home/jxs794/PROTAC_BUILDER/backup_receiver` returns:
    - `targets 20`
    - `protein 20`
    - `summary 2509`
    - `chains 8935`

## Known Issues
- I did not redeploy Heroku from this workspace in this turn, so the final browser-side loading overlay and updated metric labels still need one production smoke test after deploy.
- I validated the template JavaScript syntactically, but I did not run an in-app browser visual pass from this environment.
- One example structure-detail request I tried with a guessed virus/protein combination returned `404`; that looked like a test-input mismatch rather than a routing failure, and the main structure-detail route itself was not changed in this task.

## Manual Verification
1. Deploy the updated V-LiSEMOD app to Heroku.
2. Open `/protacability_page`.
3. Confirm the loading overlay appears immediately on page load and disappears after filters + first results finish loading.
4. Check the second summary card label in each view:
   - `Target Browser` -> `Mapped/exposed target groups`
   - `Protein Summary` -> `Mapped/exposed protein groups`
   - `Structure Summary` -> `Mapped/exposed structures`
   - `Chain Details` -> `Mapped/exposed chain rows`
5. Click the four view buttons and confirm:
   - the active button changes
   - the results layout changes
   - the summary cards update
   - pagination resets to page 1
   - the overlay appears while loading
6. Click `Apply Filters`, `Reset`, `Previous`, `Next`, and change page size; confirm the overlay appears and then clears.
7. Open a detail panel and confirm:
   - modal shows loading feedback
   - content loads
   - overlay does not get stuck on failure
8. Optional route checks:
   - `curl -sS "https://vlisemod-0e358c20a94d.herokuapp.com/api/protacability/filter_options" -o /tmp/heroku_filter_options.json`
   - `curl -sS "https://vlisemod-0e358c20a94d.herokuapp.com/api/protacability/search?view=targets&page=1&page_size=50" -o /tmp/heroku_targets.json`
   - `curl -sS "https://vlisemod-0e358c20a94d.herokuapp.com/api/protacability/search?view=chains&page=1&page_size=50" -o /tmp/heroku_chains.json`
9. Watch Heroku logs and confirm there are no new:
   - `H12 Request timeout`
   - `sqlite3.OperationalError`
   - stuck `/protacability/source` requests for normal UI loads

## Suggested Next Prompt
If the production smoke test still feels slow after deploy, use:

“Optimize the PROTACability detail and export flows so modal/detail loads and CSV exports keep response times predictable under broad filters, while preserving the current grouped summary semantics.”
