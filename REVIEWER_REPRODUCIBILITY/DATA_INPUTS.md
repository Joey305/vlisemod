# Frozen inputs and provenance

## Included small inputs

| Input | Purpose | Release content |
| --- | --- | --- |
| `manifests/FROZEN_CIF_CORPUS_MANIFEST.csv` | Corpus identity/checksum manifest | 11,533 mmCIF records; paths are relative only |
| `inputs/chemistry/frozen_component_chemistry.csv` | Explicit component chemistry for Stage 06 | 2,830 resolved component identities |
| `inputs/mapping/frozen_mapping_remediation_registry.csv` | Explicit pending mapping registry for Stage 06a | 719 occurrence-level rows |
| `fixture/PDB_FILES/.../3EKY.cif` | Small 3EKY/DR7 end-to-end fixture | two contextual corpus paths |

`INPUT_SHA256.txt` is the machine-verifiable checksum list for every included small frozen input.

The chemistry file was derived author-side from the final validated release `ligands` table and is the only chemistry source used at reviewer runtime. It replaces historical fallback reads from the parent website database. The remediation registry was derived author-side from the final release pending queue using stable occurrence keys (entry/model/label asymmetry/component/auth chain/residue/insertion code). Its expected queue distribution is 561 bond-template, 99 exact-graph, 30 local-recovery, 22 complex-graph, 5 fragment-review, and 2 curation-review rows.

## Public structural corpus

The full corpus is not placed in the ZIP and no 15 GB archive needs to be distributed by this project. The structures are public RCSB PDB mmCIF files. The package’s frozen manifest records the 11,533 exact virus/protein hierarchy paths analyzed in the release, representing 7,610 unique PDB accessions. The downloader retrieves each unique accession once, validates it, and materializes the exact recorded relative paths.

Retrieve and verify the corpus before a full run:

```bash
python reproduce.py --download-inputs --source ./PDB_FILES --workers 6
python reproduce.py --verify-inputs --source ./PDB_FILES
```

This enforces exactly 11,533 expected CIF paths and SHA-256 values. It rejects missing/modified files and rejects extras for a frozen reconstruction. The manifest SHA-256 is the release’s raw coordinate-file checksum. If an upstream accession has changed since release, the downloader records `UPSTREAM_REVISION_CHANGED` with the expected and observed hashes. The optional `--allow-current-upstream` path remains visibly distinct from frozen-release reproduction.

## Data scope

The package contains scripts and derived frozen metadata under MIT terms. Coordinate-file redistribution and any component-dictionary terms remain governed by their upstream providers; provide the external corpus through an appropriate archival route rather than assuming the software license covers those data.
