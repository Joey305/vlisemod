# Reproduction protocol

## Environment

Use the locked Conda recipe:

```bash
conda env create -f environment.yml
conda activate vlisemod-reviewer-reproducibility
python reproduce.py --check-environment
```

The check requires Python modules `gemmi`, Biopython, NumPy, SciPy, RDKit, and pandas, plus the `pdbe-arpeggio` executable and its OpenBabel dependency.

## Fast end-to-end fixture

```bash
python reproduce.py --fixture --workers 1
```

This is a real calculation against frozen 3EKY CIF input, not a mocked database. Two identical coordinate copies are kept under their distinct corpus paths so the fixture also exercises contextual classification. Stage 15 has full-corpus denominators and is therefore not run for the compact fixture.

## Public structural-input acquisition and full scientific reconstruction

Retrieve the public structural inputs from RCSB PDB, then verify them locally before running the full reconstruction:

```bash
python reproduce.py --download-inputs --source ./PDB_FILES --workers 6
python reproduce.py --verify-inputs --source ./PDB_FILES
python reproduce.py --full --source ./PDB_FILES \
  --database outputs/viral_data_cif_v2_REPRODUCED_RELEASE.db --workers 8
```

The manifest has 11,533 expected hierarchy rows and 7,610 unique PDB entries. The downloader caches each unique public entry under `.download_cache`, validates Gemmi parseability and `_entry.id`, then hard-links or copies it to each exact manifest relative path. It retries network failures, writes atomically through `.part` files, and records JSON/CSV/Markdown reports. Strict verification accepts only SHA-256 bytes matching the frozen release manifest. It never downloads during `--verify-inputs`.

If RCSB has revised an entry, normal verification fails deliberately and the report marks it `UPSTREAM_REVISION_CHANGED`. `--allow-current-upstream` permits a parsed but checksum-different hierarchy and labels a later run `CURRENT_UPSTREAM_RECONSTRUCTION`; it does not relax Stage 15 expectations.

The orchestration verifies all 11,533 manifest checksums before creating a database. It executes Stages 01, 02, 03, 05, 06, 06a, 07–15. Stage 07 retries deterministic timeout rows at a longer timeout; Stage 09 retries incomplete included cases with `--resume` and fails when a retry makes no progress. Stage 15 remains strict: any failed or warning check is release-blocking according to the final validator.

After Stage 15, the orchestrator writes `outputs/SCIENTIFIC_CONTENT_FINGERPRINTS.json` and compares it to the frozen reference. `outputs/REPRODUCED_RELEASE_SHA256.txt` records the produced database, validation reports, and fingerprints.

## Partial and recovery runs

`--from-stage` and `--to-stage` select an inclusive canonical range; `--resume` passes resume semantics to Stages 07 and 09.

```bash
python reproduce.py --full --source /path/to/PDB_FILES --database outputs/recovery.db \
  --from-stage mapping --to-stage validate --resume --workers 8
```

Use the same database and frozen inputs for recovery. A partial run that does not finish Stage 15 is not a releasable reconstruction.

`--dry-run` prints the exact commands and never creates a database. `--limit` is intended only for bounded development diagnostics; it cannot satisfy Stage 15.

## Optional website compatibility augmentation

Only after a successful scientific full run:

```bash
python reproduce.py --web-compat --database outputs/viral_data_cif_v2_REPRODUCED_RELEASE.db
```

This copies the scientific database to `outputs/viral_data_cif_v2_WEB_COMPAT.db` before running website-facing synonym/diagram builders. The scientific input is never modified.
