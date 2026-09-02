# Package audit record

## Build checks

- Canonical scientific scripts are present for Stages 01, 02, 03, 05, 06, 06a, and 07–15.
- Historical Stage 04 wrapper and old geometry/attachment scripts are absent from the canonical stage list.
- Reviewer runtime source has no dependency on a parent website database or parent chemistry file.
- The package manifest contains only relative corpus paths; all included small inputs and canonical scripts have SHA-256 manifests.
- The distributed ZIP excludes `outputs/`, databases, logs, Python bytecode, and private/deployment markers.

## Executed checks

- Unit tests: 10/10 passing (`tests/test_reproduce.py`, `tests/test_download_inputs.py`).
- Environment check: Python 3.13.12; required modules and `pdbe-arpeggio` found.
- Final database release validation (read-only against frozen final release): PASS, 0 failed checks, 0 warnings.
- Scientific-content fingerprint comparison against the frozen final release: PASS.
- ZIP extraction audit: PASS.
- Extracted-ZIP 3EKY/DR7 end-to-end fixture: PASS through real Stages 01–13.
- Live RCSB 3EKY downloader, resume, and local verification: PASS; current RCSB bytes match the frozen manifest and local release corpus.

The ZIP contains no general-release CIF files. Its only coordinate files are the two explicit fixture copies: `fixture/PDB_FILES/HIV_1/capsid_protein/3EKY.cif` and `fixture/PDB_FILES/HIV_1/protease/3EKY.cif`.

The full 11,533-CIF reconstruction was intentionally not rerun during package qualification. Public RCSB structures are retrieved on demand from the frozen manifest; the exact full command and its strict Stage-15/fingerprint checks are documented in `REPRODUCIBILITY.md`.
