# Qodex.summary

## Task
Migrate V-LiSEMOD ligand image and PyMOL session database reads to RANDY.

## Original Goal
Fix Heroku failures where `/generate_ligand_images` and `/generate_pymol_session` still tried to query local SQLite tables that do not exist on Heroku.

## Assumptions
- `VLISMOD_DATA_BACKEND=randy` is the intended Heroku runtime mode for database-backed V-LiSEMOD routes.
- The production RANDY base is `https://randy.rove-vernier.ts.net/backup/vlismod`.
- RANDY reads the production database from `/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`.
- The Heroku dyno should still generate SVG and `.pml` files locally, but should not perform database reads locally.
- `backup_receiver.app:APP` is the live RANDY service that must expose the new V-LiSEMOD POST endpoints.

## Files Inspected
- `app.py`
- `RANDY/vlismod_data_routes.py`
- `RANDY/app.py`
- Remote: `/home/jxs794/PROTAC_BUILDER/backup_receiver/app.py`
- Remote: `/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`
- Heroku logs provided in the prompt

## Files Changed
- `app.py`
- `RANDY/vlismod_data_routes.py`
- `Qodex.summary.md`
- Remote runtime sync for validation: `/home/jxs794/PROTAC_BUILDER/backup_receiver/vlismod_data_routes.py`

## Files Created
- None

## Implementation Summary
Two new authenticated RANDY endpoints were added under both `/api/vlismod` and `/backup/vlismod`:
- `POST /ligand-images-data`
- `POST /pymol-session-data`

`/generate_ligand_images` now fetches its SMILES rows, chain/residue rows, and solvent-exposed SMILES atom indices from RANDY in strict `randy` mode, then continues generating SVG files locally.

`/generate_pymol_session` now fetches ligand-chain resolution, functional-group atoms, binding-pocket rows, distal atoms, solvent-exposed atoms, hydrated atoms, and Rupley SASA atoms from RANDY in strict `randy` mode, then continues writing the `.pml` file locally.

Local fallback behavior was preserved for `local` and `auto` modes, but strict `randy` mode now returns clear remote-failure errors instead of silently touching local SQLite.

## Key Decisions
- Heroku should not read local SQLite because the dyno does not have the production V-LiSEMOD database.
- SVG generation and PyMOL script generation can remain local because they are file-generation tasks, not database-access tasks.
- Database reads were moved to RANDY because RANDY already hosts the real production database and authenticated backup URL.
- Strict `randy` mode now preserves meaningful remote statuses such as `404`, and returns `502` only for actual RANDY request failures.

## Commands Run
- Code inspection with `grep`, `sed`, and `rg`
- Local schema/value checks with `sqlite3`
- `python -m py_compile app.py RANDY/app.py RANDY/vlismod_data_routes.py`
- Local Flask test-client checks for RANDY blueprint routes
- Local strict `randy`-mode app checks with mocked RANDY POST responses
- Local strict `randy`-mode app checks against the live RANDY HTTPS endpoint
- Remote RANDY runtime inspection with `ssh`
- Remote file sync of `vlismod_data_routes.py` into `backup_receiver`
- Remote Gunicorn worker reload with `kill -HUP <master-pid>`
- Live HTTPS `curl` checks against `https://randy.rove-vernier.ts.net/backup/vlismod`
- Git safety checks with `git status --short` and `rg`

## Validation Results
- Root cause confirmed:
  - `/generate_ligand_images` still queried `Functional_GROUPED`, `Ligand_Arp_Diagram`, `RUPLEY_SASA_DATA`, and `SMILES_MAP_PDB` locally.
  - `/generate_pymol_session` still queried `ligand_atoms`, `Functional_Group_Atoms`, `receptor_binding_pocket`, `distal_atoms`, and `RUPLEY_SASA_DATA` locally.
  - This matched the Heroku `sqlite3.OperationalError: no such table ...` failures from the prompt.

- Syntax:
  - `python -m py_compile app.py RANDY/app.py RANDY/vlismod_data_routes.py` passed.

- RANDY endpoint checks:
  - Internal `http://127.0.0.1:8787/backup/vlismod/health` -> `200`
  - External `https://randy.rove-vernier.ts.net/backup/vlismod/db-health` -> `200`
  - External unauthenticated `POST /backup/vlismod/pymol-session-data` -> `401`
  - External authenticated `POST /backup/vlismod/ligand-images-data` for `Human immunodeficiency virus 1 / 3EKY / DR7 / A` -> `200`
  - External authenticated `POST /backup/vlismod/pymol-session-data` for `3EKY / DR7 / A` -> `200`
  - `db-health` confirms:
    - `db_exists=true`
    - `db_path=/home/jxs794/PROTAC_BUILDER/VLISEMOD/Database/viral_data.db`
  - `ligand-images-data` returned:
    - `smiles_data` count `1`
    - `chain_residue_data` count `1`
    - `solvent_exposed_atom_map` key `3EKY|DR7|A`
  - `pymol-session-data` returned:
    - `ligand_chain=A`
    - `functional_groups` count `16`
    - `binding_pocket` count `121`
    - `distal_atoms` count `6`
    - `solvent_exposed_atoms` count `12`
    - `hydrated_atoms` count `51`
    - `rupley_sasa` count `12`

- V-LiSEMOD route checks:
  - Strict `randy` mode with mocked RANDY data:
    - `/generate_ligand_images` -> `200`
    - `/generate_pymol_session` -> `200`
  - Strict `randy` mode with mocked RANDY failure:
    - `/generate_ligand_images` -> `502`
    - `/generate_pymol_session` -> `502`
  - Strict `randy` mode against live RANDY HTTPS:
    - `/generate_ligand_images` for `Human immunodeficiency virus 1 / 3EKY / DR7 / A` -> `200`
    - `/generate_pymol_session` for `3EKY / DR7 / A` -> `200`

- Heroku log verification:
  - Existing logs in the prompt correctly showed the pre-fix local SQLite failure mode.
  - A fresh Heroku log verification for these two routes was not possible yet because this turn did not push or redeploy Heroku.

## Known Issues
- Heroku itself has not been redeployed with these new code changes yet, so the live Heroku app may still show the old failures until the next deploy.
- Many other complex V-LiSEMOD routes still query local SQLite and may fail later in strict `randy` mode until they are migrated too.
- The live RANDY runtime file was updated directly for validation, so the local repo and remote runtime are now aligned for `vlismod_data_routes.py`, but Heroku still needs the matching repo deploy.

## Manual Verification
After the next V-LiSEMOD deploy, verify:

```bash
curl -i https://vlisemod-0e358c20a94d.herokuapp.com/get_viruses
curl -i https://vlisemod-0e358c20a94d.herokuapp.com/get_ligands/3EKY
curl -i https://vlisemod-0e358c20a94d.herokuapp.com/check_functional_groups/3EKY
```

In the browser:
1. Open the homepage workflow.
2. Select `Human immunodeficiency virus 1`.
3. Select PDB `3EKY`.
4. Select ligand `DR7`.
5. Generate ligand images.
6. Generate the PyMOL session.

Expected:
- ligand images page renders successfully
- the solvent-exposed SVG is generated
- the PyMOL `.pml` download succeeds
- Heroku logs no longer show `no such table: Functional_GROUPED` or `no such table: ligand_atoms` for those requests

Direct RANDY checks:

```bash
TOKEN="<token from env only>"
BASE="https://randy.rove-vernier.ts.net/backup/vlismod"

curl -i -H "Authorization: Bearer $TOKEN" "$BASE/db-health"
curl -i -X POST \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"virus_name":"Human immunodeficiency virus 1","pdb_code":"3EKY","ligand_name":"DR7","chain":"A"}' \
  "$BASE/ligand-images-data"

curl -i -X POST \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"pdb_code":"3EKY","ligand_name":"DR7","chain":"A","options":{"functional_groups":true,"binding_pocket":true,"distal_atoms":true,"solvent_exposed_atoms":true,"hydrated_atoms":true,"rupley_sasa":true}}' \
  "$BASE/pymol-session-data"
```

## Suggested Next Prompt
Migrate the next batch of strict-`randy` V-LiSEMOD routes that still use local SQLite, especially chart generation, ligand interaction comparison, SASA-heavy visualizations, and PROTACability endpoints, while preserving local file generation and keeping Heroku stateless for database reads.
