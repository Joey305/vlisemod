# V-LiSEMOD CIF-only database rebuild

This directory is a self-contained, non-production rebuild pipeline. It reads only the retained `PDB_FILES` mmCIF corpus. It never downloads from RCSB, never edits CIFs, and never writes to `viral_data.db`.

## Foundation model

`structures` represents an exact retained CIF revision (`entry_id + SHA-256`). Folder labels are stored separately in `structure_classifications`. `ligands` is reusable chemistry; `ligand_instances` is a physical occurrence distinguished by the full mmCIF tuple: structure revision, deposited model number, label asym/component, author asym/residue and normalized insertion code. `ligand_instance_atoms` retains all deposited atom-site rows, including alternate locations.

Insertion `.` and `?` normalize to empty string only for comparisons; their raw values remain stored. IDs are all SQLite `TEXT` with no fixed-length assumptions. The cohesive conformer policy selects blank/shared atoms plus one named altloc by heavy-atom completeness, summed occupancy, then lexical altloc ID.

## Setup

Use an isolated environment and install the exact project requirements:

```bash
cd /Users/jxs794/Documents/VLISEMOD/CIF_DATABASE_REBUILD
python -m pip install -r requirements.txt
```

`pdbe-arpeggio` is deliberately not listed: it requires separate installation and selector validation before the Arpeggio stage is enabled.

## Arpeggio recovery and provenance

The Arpeggio stage uses a bounded, persisted retry policy. A normal first pass
uses the authoritative source mmCIF (or the existing canonical single-model
derived input when model/altloc handling requires it). `--resume` always skips
any ligand instance having a valid completed run.

For a recorded parser, duplicate-identity, or OpenBabel/BioPython mapping
failure, the next resume attempt uses a deterministic `sanitized_full` mmCIF.
For a recorded timeout, it uses a `sanitized_pocket` input containing complete
residues having any atom within 12 Angstrom of the selected ligand. This radius
is deliberately larger than Arpeggio's contact-search neighborhood. A timeout
retry defaults to 600 seconds and is attempted once per resume invocation; a
subsequent failure remains explicit and finite.

Derived inputs never replace or edit source mmCIF files. Each derived attempt
is stored below `outputs/arpeggio/raw/<run_id>/<strategy>/` with:

- `derived_input.cif`, using one deposited model, coherent residue-wide altloc
  selection, sequential derived atom IDs, and safe CIF quoting;
- `atom_provenance.jsonl`, mapping every retained derived atom to its original
  atom-site identity and the canonical ligand atom where applicable;
- `derivation_manifest.json`, recording source checksum, selector,
  deterministic operations, context radius, and derived checksum;
- stdout, stderr, and validated Arpeggio JSON output.

Attempt details are persisted in `arpeggio_attempts`, retained-atom mappings in
`arpeggio_derived_atom_map`, and final per-instance summaries in
`ligand_arpeggio_runs`. A subprocess exit code of zero is insufficient by
itself: JSON must be present and parseable, the selected ligand must be present,
and ligand endpoints must reconcile to canonical atoms (or to an explicitly
recognized multi-atom aromatic group).

Classify the original unresolved set without rerunning it:

```bash
python 29_arpeggio_failure_audit.py --database ./viral_data_cif_v2.db --baseline-only
```

Resume unresolved instances only:

```bash
python run_pipeline.py --stage arpeggio --database ./viral_data_cif_v2.db \
  --resume --workers 10 --per-instance-timeout 300 --retry-timeout 600 \
  --fallback-radius 12 --progress-every 25
```

The first recovery invocation sends prior parser failures to `sanitized_full`
and prior timeouts to `sanitized_pocket`. If a sanitized full attempt itself
times out, the persisted status causes the next `--resume` invocation to use
the controlled pocket retry; successful prior results are never regenerated.

## Commands

Create a fresh development database (remove only this development DB first if a true fresh rebuild is wanted):

```bash
python 02_create_database.py --database ./viral_data_cif_v2.db
```

Freeze the retained corpus:

```bash
python run_pipeline.py --stage inventory
```

Ingest all structures, classifications, candidates, and atom sites:

```bash
python run_pipeline.py --stage ingest --resume
```

Run the reviewed foundation sequence:

```bash
python run_pipeline.py --all --resume
```

Restrict any foundation command with `--pdb-id 5Y9E` or `--limit 10`. `--workers` is accepted for the future expensive stages; foundation ingestion is deliberately transactional/single-writer. `--dry-run` shows selected stages. Re-runs use unique constraints and upserts instead of destructive replacement.

Unit/fixture test:

```bash
python -m unittest discover -s tests -v
```

## Outputs and milestones

Inventory outputs live in `manifests/`; foundation results are in `viral_data_cif_v2.db` and `outputs/FOUNDATION_VALIDATION_REPORT.md`. Mapping, SASA, functional-group, geometry, Arpeggio, PROTACability, attachment, compatibility and legacy-comparison entry points are present but intentionally stop until their specified validation milestones are approved. The preserved legacy method boundaries are documented in the parent audit: CIF-native input/identity changes are required; scientific settings such as Shrake–Rupley (1.40 Å, atom-level, >0.1 Å² derived exposure), RDKit ETKDG/UFF/MCS settings, existing SMARTS, Arpeggio contact engine, and PROTACability/attachment scoring are to remain unchanged for validation comparison.
