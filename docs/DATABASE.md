# Database Guide

## Core Database Purpose

The main application database is `viral_data.db`, a SQLite store used to support structure lookup, ligand-centered annotations, interaction analysis, export workflows, and PROTACability triage. The app reads directly from this database in many routes, and some enrichment/import scripts regenerate or append data layers outside the request cycle.

Unless otherwise noted, table names below were confirmed from the repository’s current `viral_data.db`. When a workflow clearly expects a table but the surrounding conventions are inconsistent, the table is described as used or expected by the workflow.

## Viral Protein Metadata

- Confirmed table: `Virus_Proteins`
- Purpose: virus names, protein labels/types, and PDB-linked target metadata that drive target-centric query workflows.
- Route usage: used by the protein-query workflows and distinct-option loaders such as `/get_virus_names_list_distinct`, `/get_protein_types_list_distinct`, and `/get_pdbs_for_virus_protein`.

## Ligand Synonyms and Identifiers

- Confirmed table: `ligand_synonyms`
- Route usage: ligand lookup endpoints query `Ligand_Synonyms` in code, which SQLite resolves case-insensitively against `ligand_synonyms`.
- Purpose: map ligand 3-letter codes to synonyms for ligand-first selection, comparison, and query filtering.

## Ligand Atoms

- Confirmed table: `ligand_atoms`
- Purpose: atom-level coordinates and atom identities for ligand instances in viral structures.
- Typical fields: `virus_name`, `pdb_id`, `ligand`, `chain`, `atom_id`, `exact_atom`, `atom_type`, `x`, `y`, `z`.
- Workflow role: structure explorer and export surfaces use this as a core atom-level ligand layer.

## Ligand SMILES / Atom Mappings

- Confirmed tables: `Ligand_Atoms_Smiles`, `SMILES_MAP_PDB`
- Purpose:
  - `Ligand_Atoms_Smiles` stores ligand-level SMILES and derived properties such as molecular weight and functional-group summaries.
  - `SMILES_MAP_PDB` stores mapping between PDB atom identities and SMILES atom identities.
- Workflow role: ligand comparison, SVG rendering, mapped-atom logic, and warhead/linkability scoring.

## Arpeggio Contacts

- Confirmed tables: `Arpeggio_Contacts_Data`, `Ligand_Arp_Diagram`
- Purpose:
  - `Arpeggio_Contacts_Data` stores atom/residue contact calls and distances.
  - `Ligand_Arp_Diagram` stores ligand/PDB/chain/ligand-instance linkage used by interaction display workflows.
- Workflow role: interaction charts, ligand comparison, and structure interpretation.

## Ligand Diagrams

- Confirmed table: `Ligand_Arp_Diagram`
- Generated/static companions: the app also generates chart images and ligand SVG-like outputs in static/output folders at runtime.
- Workflow role: 2D ligand interaction or chart-driven display workflows.

## SASA / Solvent Exposure Data

- Confirmed tables: `RUPLEY_SASA_DATA`, `solvent_exposed_atoms`
- Related table: `ligand_water_distances`
- Purpose:
  - `RUPLEY_SASA_DATA` stores Shrake-Rupley solvent accessibility results used in the main ligand view and PROTACability-linked atom interpretation.
  - `solvent_exposed_atoms` appears to store a summarized atom-level exposed-atom layer used in exports/workflows.
  - `ligand_water_distances` appears to support solvent/water-context analysis where present.
- Workflow role: surface exposure review and modifiable-atom reasoning.

## Functional Group Annotations

- Confirmed tables: `Functional_GROUPED`, `Functional_Group_Atoms`
- Purpose:
  - `Functional_GROUPED` stores ligand-level grouped functional annotations.
  - `Functional_Group_Atoms` stores atom-level functional-group membership.
- Workflow role: optional annotations in the structure explorer and a source layer for warhead/linkability interpretation.

## Binding Pocket Data

- Confirmed table: `receptor_binding_pocket`
- Purpose: protein atoms or residues within the defined ligand-adjacent binding context.
- Workflow role: optional binding-pocket object generation and query/export support.

## Additional Ligand Context Layers

- Confirmed tables: `distal_atoms`, `Covalent_Noncovalent`
- Purpose:
  - `distal_atoms` appears to capture ligand atoms distal from nearby protein atoms under the project’s geometric workflow.
  - `Covalent_Noncovalent` classifies inhibitor type at the ligand-structure level.
- Workflow role: structural interpretation and enrichment rather than primary navigation.

## PROTACability Tables

- Confirmed tables:
  - `protacability_assessment`
  - `protacability_lysine_proximity`
  - `protacability_ligand_inventory`
  - `protacability_warhead_linkability`
  - `protacability_degrader_readiness`
- Purpose:
  - `protacability_assessment`: summary structural-priority table used by the dashboard.
  - `protacability_lysine_proximity`: lysine accessibility and ligand-proximity cues.
  - `protacability_ligand_inventory`: ligand instances present in each structure.
  - `protacability_warhead_linkability`: ligand-centered linkerability and atom-attachment evidence.
  - `protacability_degrader_readiness`: combined readiness heuristics and geometry cues.
- Import workflow: `TOOLS/import_protacability_data.py` imports CSV outputs from `PDB_FILES/` into these tables and adds indices.

## Generated / Static Data Assets

The repo references or generates several filesystem-backed assets outside the database:

- `PDB_FILES/`: expected location for imported/generated PROTACability CSVs and likely coordinate-related source material
- `output_csvs/`: temporary export products
- `pml_sessions/`: generated PyMOL-oriented session artifacts
- `static/charts/`: generated comparison or interaction charts
- `static/coordinate_cache/`: cached coordinate files for viewers
- `static/ligand_sdf_cache/`: cached ligand instance SDF files
- `static/ligand_images/`: generated ligand imagery

These assets are operationally important but should generally be treated as generated runtime artifacts rather than canonical source data.

## Data That Should Not Be Committed Publicly

Do not commit the following to public GitHub unless a deliberate redacted sample is created:

- production or local copies of `viral_data.db`
- `users.db` and any user/session records
- `PDB_FILES/` bulk structure or generated CSV bundles
- generated caches in `static/coordinate_cache/`, `static/ligand_sdf_cache/`, `static/charts/`, and ligand-image folders
- PyMOL sessions, exports, and temporary outputs
- local model weights, checkpoints, and any token-bearing configuration

## Schema Confidence Notes

- Table presence was verified from the current local SQLite database.
- Route usage was inferred from the current Flask codebase.
- Some naming is historically inconsistent in code (`Ligand_Synonyms` vs `ligand_synonyms`), but SQLite’s case-insensitive table resolution makes the current workflow function.
