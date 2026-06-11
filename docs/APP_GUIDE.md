# Application Guide

## Overview

V-LiSEMOD is organized around structure-guided viral ligand analysis. Users can begin from a structure, a target, or a ligand and move toward solvent exposure review, interaction interpretation, and degrader-readiness triage without account creation or login in the current public-facing workflow.

## Home / Structure Explorer

- Purpose: main entry point for selecting a virus, PDB structure, and ligand, then generating PyMOL-oriented outputs and ligand views.
- Typical user question: "What does this ligand look like in a specific viral co-crystal structure, and which atoms appear exposed enough to modify?"
- Main capabilities: virus/PDB/ligand selection, chain-aware selection, PyMOL session generation, ligand image generation, optional binding-pocket export, SASA highlight option, and functional-group visualization when available.
- Major data dependencies: `ligand_atoms`, `Functional_Group_Atoms`, `receptor_binding_pocket`, `distal_atoms`, `RUPLEY_SASA_DATA`, `Ligand_Arp_Diagram`, `Functional_GROUPED`.
- Public interpretation note: describe this as a structure exploration and export workspace, not as an automated design engine.

## About

- Purpose: scientific framing, workflow explanation, and companion-tool context.
- Typical user question: "Why does V-LiSEMOD exist and how should it be used?"
- Main capabilities: project rationale, module framing, and ecosystem links.
- Major data dependencies: template-driven page with limited direct data dependency.
- Public interpretation note: this page should reinforce that V-LiSEMOD is evidence-oriented and hypothesis-generating.

## Protein Query

- Purpose: target-centric filtering and export assembly.
- Typical user question: "Which structures exist for this viral target, and can I gather the linked data layers for review?"
- Main capabilities: filter by virus name, protein type, optional ligand, matched PDB selection, export-ready dataset selection, and PROTACability dashboard handoff.
- Major data dependencies: `Virus_Proteins`, `ligand_synonyms`, and the export-linked source tables used by the query workflow, including optional PROTACability tables when imported.
- Public interpretation note: frame this module as a curated query and export surface, not as an arbitrary schema browser.

## Ligand Indexer

- Purpose: ligand-first lookup across mapped structural contexts.
- Typical user question: "Where does this ligand appear, and what interaction chart can I generate for a chosen instance?"
- Main capabilities: ligand or synonym search, mapped PDB/chain/residue selection, interaction chart generation, and interaction-type definitions.
- Major data dependencies: `ligand_synonyms`, `Ligand_Arp_Diagram`, `Arpeggio_Contacts_Data`.
- Public interpretation note: this is a practical lookup and charting workflow, not an exhaustive ligand-similarity engine.

## Ligand Comparison

- Purpose: compare one ligand across multiple structural contexts.
- Typical user question: "How does this ligand’s interaction behavior change across different viral protein structures?"
- Main capabilities: multi-context ligand selection, comparison charts, ligand SVG display, atom-level interaction burden review, and external handoff to PROTAC Builder.
- Major data dependencies: `ligand_synonyms`, `Ligand_Atoms_Smiles`, `SMILES_MAP_PDB`, `Arpeggio_Contacts_Data`.
- Public interpretation note: useful for cross-structure comparison and exit-vector discussion, but still limited by structure quality and mapping coverage.

## PROTACability Assessment

- Purpose: transparent structural-priority and degrader-readiness triage.
- Typical user question: "Which viral target structures look more promising for degrader-oriented follow-up, and why?"
- Main capabilities: target, protein, structure, and chain-detail views; filters for tiers and evidence layers; CSV export; coordinate and ligand-instance support for structural inspection.
- Major data dependencies: `protacability_assessment`, `protacability_lysine_proximity`, `protacability_ligand_inventory`, and when present `protacability_warhead_linkability` plus `protacability_degrader_readiness`.
- Public interpretation note: results are heuristic, transparent, and intended for triage. They are not experimentally validated degradation predictions.

## PROTAC Builder Handoff

- Purpose: send promising ligand contexts into a downstream degrader-design environment.
- Typical user question: "I found a plausible warhead context; where do I continue linker and recruiter design work?"
- Main capabilities: external links and legacy redirect compatibility to a configurable PROTAC Builder URL.
- Major data dependencies: no local analytic dependency beyond the selected ligand context and configured external URL.
- Public interpretation note: PROTAC Builder is a companion tool, not a module maintained inside this repository.

## Optional Drug GPT / BioGPT Module

- Purpose: optional embedded assistant workflow when feature flags enable the `/drugapp/` blueprint.
- Typical user question: "Can I ask a local assistant about a molecule or context from inside the app?"
- Main capabilities: assistant page and query route when `ENABLE_DRUG_GPT=1`; disabled-state placeholder routes otherwise.
- Major data dependencies: local model configuration, Hugging Face token support, and runtime packages imported by `DRUGapp.py`.
- Public interpretation note: describe this module as optional and deployment-dependent. In the default lightweight public mode it should be considered disabled unless explicitly provisioned.
