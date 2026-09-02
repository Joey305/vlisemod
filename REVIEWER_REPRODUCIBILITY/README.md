# V-LiSEMOD reviewer reproducibility package

This package reproduces the frozen CIF-native scientific database release. It is deliberately separate from the Flask website and does not read, alter, or require a website database.

## Quick start

```bash
conda env create -f environment.yml
conda activate vlisemod-reviewer-reproducibility
python reproduce.py --check-environment
python reproduce.py --fixture --workers 1
```

The fixture runs the real current Stages 01–13 on the two frozen contextual copies of `3EKY.cif`, then verifies DR7 at chain A/residue 100: 51 mapped atoms, 12 solvent-exposed atoms, target scores A=52.86 and B=81.43, no High site, and conditional Moderate atoms CAO/CAR/CAS/NBD (score 79).

## Three review levels

1. **Fast fixture (minutes):** `python reproduce.py --fixture`.
2. **Retrieve and verify public structural inputs:** `python reproduce.py --download-inputs --source ./PDB_FILES`, then `python reproduce.py --verify-inputs --source ./PDB_FILES`.
3. **Full scientific reconstruction:** `python reproduce.py --full --source /path/to/PDB_FILES --database outputs/viral_data_cif_v2_REPRODUCED_RELEASE.db --workers 8`.

The package retrieves public mmCIF structures directly from RCSB PDB. Its frozen manifest records 11,533 classification paths representing 7,610 unique PDB entries; repeated accessions are downloaded once then materialized at each exact manifest path. Current RCSB byte revisions are detected rather than silently accepted. The full command runs the canonical stages through strict Stage 15 and compares content fingerprints with the release reference. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) and [DATA_INPUTS.md](DATA_INPUTS.md).

## Scope boundary

`--full` is the manuscript scientific pipeline. It uses only the supplied frozen chemistry and mapping-registry inputs plus the frozen CIF corpus. `--web-compat` is a separate, optional post-release augmentation that writes a copy with website-facing synonym/diagram tables; it is not needed for any scientific count, score, or Stage-15 check.

## Canonical stages

`01 → 02 → 03 → 05 → 06 → 06a → 07 → 08 → 09 → 10 → 11 → 12 → 13 → 14 → 15`

Stage 04 is intentionally absent because it is a historical compatibility wrapper around ingestion. Old geometry and attachment scripts are excluded. Full-release versions are documented in [METHODS_AND_VERSIONS.md](METHODS_AND_VERSIONS.md).

## Integrity files

- `INPUT_SHA256.txt` covers frozen small inputs and fixture CIF files.
- `SCRIPT_SHA256.txt` covers the orchestrator and canonical stages.
- `reference/SCIENTIFIC_CONTENT_FINGERPRINTS.json` covers stable scientific content from the final release, excluding IDs, timestamps, and local paths.

The distributed ZIP excludes `outputs/`; all run databases, logs, and derived reports are generated locally.
