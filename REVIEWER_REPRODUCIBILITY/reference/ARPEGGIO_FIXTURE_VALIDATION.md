# Arpeggio fixture validation

- Runtime: `<author-home>/miniconda3/envs/viraldb2/bin/pdbe-arpeggio` 1.4.4.
- Input provenance: each run records the source checksum; derived analysis-only
  inputs additionally record their path, checksum, and canonical deposited
  model. Source CIF files are never modified.
- Endpoint reconciliation is occurrence-scoped by `ligand_instance_id` and
  canonical selected conformer. Raw interaction labels and deduplicated
  ligand/environment atom pairs remain separate database records.

| Fixture | Ligand instance | Canonical identity | Input / selector | Raw labels | Unique pairs | Selected endpoints (observed/reconciled/unmatched/ambiguous) | Sibling leakage |
| --- | ---: | --- | --- | ---: | ---: | --- | ---: |
| 3EKY | 59418 | DR7; auth A/100, label C | source; `/A/100/` | 1364 | 455 | 556 / 556 / 0 / 0 | 0 |
| 5Y9E | 454181 | GOL; auth A/501, label F | source; `/A/501/` | 181 | 59 | 59 / 59 / 0 / 0 | 0 |
| 1AI1 | 331025 | AIB; auth P/323, label C | source; `/P/323/` | 295 | 48 | 48 / 48 / 0 / 0 | 0 |
| 2M3L | 453494 | ZN; canonical model 10, auth A/201, label B | derived, model-isolated analysis CIF; `/A/201/` | 221 | 27 | 27 / 27 / 0 / 0 | 0 |
| 6MCF | 13918 | RY; auth/label A/71 | source; `/A/71/` | 1020 | 306 | 306 / 306 / 0 / 0 | 0 |
| 6MCF | 13919 | RY; auth/label A/75 | source; `/A/75/` | 708 | 243 | 243 / 243 / 0 / 0 | 0 |

## Controls

- 2M3L uses a deterministic derived input containing only canonical deposited
  model 10. It is represented as analysis model 1 solely for pdbe-arpeggio
  compatibility; the deposited model number 10 is retained in provenance.
- The derived-input policy keeps shared atoms and makes a single coherent named
  conformer choice by occupancy then lexical altloc order. It never constructs
  a per-atom hybrid conformer.
- 5Y9E confirms the CLI selector namespace is auth rather than label: its label
  chain is F but successful selector is `/A/501/`.
- 1AI1's requested occurrence has no ligand insertion code; the canonical
  identity fields are nevertheless retained in the run record. The derived
  writer preserves insertion-code fields whenever it is required.
- Both 6MCF RY siblings were analysed independently and endpoint lookup is
  restricted to the requested occurrence's selected atom records.

ARPEGGIO FIXTURES PASSED — BULK RUN APPROVED
