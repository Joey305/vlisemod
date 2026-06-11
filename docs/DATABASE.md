# Database Guide

## Core Database Purpose

The main application database is `viral_data.db`, a SQLite store used to support structure lookup, ligand-centered annotations, interaction analysis, export workflows, and PROTACability triage. The app reads directly from this database in many routes, while some enrichment or import scripts regenerate data layers outside the request cycle.

The current inspected local database contains 20 tables, including the core structure, ligand, contact, SASA, and PROTACability layers described below.

## Confirmed Core Tables

### Viral Protein Metadata

- Table: `Virus_Proteins`
- Purpose: virus names, protein labels, and PDB-linked target metadata that drive target-centric workflows.

### Ligand Synonyms and Identifiers

- Table: `ligand_synonyms`
- Purpose: ligand-code and synonym mapping for ligand-first lookup, comparison, and query filtering.

### Ligand Atoms and Structural Context

- Tables: `ligand_atoms`, `distal_atoms`, `receptor_binding_pocket`
- Purpose: atom-level ligand coordinates plus ligand-adjacent structural context used by the structure explorer and export surfaces.

### Ligand SMILES and Mapping Layers

- Tables: `Ligand_Atoms_Smiles`, `SMILES_MAP_PDB`
- Purpose: ligand SMILES, derived descriptors, and atom mapping between PDB and ligand representations.

### Interaction Layers

- Tables: `Arpeggio_Contacts_Data`, `Ligand_Arp_Diagram`
- Purpose: atom/residue contact calls, interaction summaries, and ligand-structure linkage used by ligand lookup and comparison workflows.

### Surface and Solvent Exposure Layers

- Tables: `RUPLEY_SASA_DATA`, `solvent_exposed_atoms`, `ligand_water_distances`
- Purpose: solvent accessibility and related context used for exposed-atom interpretation.

### Functional Group Layers

- Tables: `Functional_GROUPED`, `Functional_Group_Atoms`
- Purpose: grouped and atom-level functional-group annotations used in structure review and warhead-linkability interpretation.

### Additional Context

- Table: `Covalent_Noncovalent`
- Purpose: ligand-structure classification context used in enrichment and interpretation workflows.

## Confirmed PROTACability Tables

- `protacability_assessment`
- `protacability_lysine_proximity`
- `protacability_ligand_inventory`
- `protacability_warhead_linkability`
- `protacability_degrader_readiness`

These tables support the PROTACability dashboard, evidence layers, detail pages, and export routes. They should be interpreted as structural triage layers rather than experimental degradation evidence.

## Import and Regeneration Notes

The repository includes `TOOLS/import_protacability_data.py` to import CSV outputs from `PDB_FILES/` into the PROTACability tables. Several other scripts in the repo also assume a local `viral_data.db` for curation, enrichment, or maintenance tasks.

## Local SQLite and RANDY Modes

V-LiSEMOD can operate in:

- local SQLite mode, where the app queries `viral_data.db` directly, or
- RANDY-backed mode, where selected route groups proxy structure and PROTACability data through the RANDY service.

In both cases, the scientific data layer is provisioned separately from the public source repository.

## Generated and Filesystem-Backed Assets

The application also references or generates several non-database assets:

- `PDB_FILES/`
- `output_csvs/`
- `pml_sessions/`
- `static/charts/`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- `static/ligand_images/`

These are operationally important but should usually be treated as generated runtime artifacts rather than canonical source material.

## Public and Private Data Boundary

Do not commit the following to public GitHub unless a deliberate redacted sample is prepared:

- local or production copies of `viral_data.db`
- legacy `users.db` files
- bulk structure or generated CSV bundles in `PDB_FILES/`
- generated caches in `static/coordinate_cache/`, `static/ligand_sdf_cache/`, `static/charts/`, and `static/ligand_images/`
- PyMOL sessions, exports, and temporary outputs
- local model weights or token-bearing configuration

## Interpretation Caution

Database presence supports present-tense claims about what the app can display or compute. It does not by itself validate the biological conclusions a user may draw from those data layers. Manuscript language should keep that distinction explicit.
