# Manuscript Outline

## Working Title

V-LiSEMOD: a structure-guided viral ligand solvent-exposed moiety database for ligand modification analysis and degrader-readiness triage

## Abstract Scaffold

Viral protein-ligand co-crystal structures contain valuable information for inhibitor optimization, exit-vector analysis, and early degrader-oriented hypothesis generation, but the relevant evidence layers are often fragmented across structural files, interaction outputs, atom annotations, and custom downstream scripts. V-LiSEMOD was developed as a web-based structural bioinformatics platform that integrates curated viral protein metadata, ligand identifiers, atom-level interaction context, solvent-exposed ligand atom analysis, functional-group annotations, and transparent PROTACability-style heuristics into a single exploration environment. The platform supports target-centric query, ligand-first indexing, cross-structure ligand comparison, PyMOL-oriented export, and degrader-readiness triage without overstating predictive certainty. V-LiSEMOD is intended as a hypothesis-generation and prioritization resource for antiviral structural biology and induced-proximity-inspired design workflows.

## Introduction Rationale

- viral structural datasets are growing, but design-relevant ligand interpretation remains fragmented
- medicinal chemistry decisions often require structure, contact, and solvent-exposure context simultaneously
- degrader-oriented thinking introduces additional target-side questions such as exposed lysine context and linker-attachment feasibility
- a transparent evidence platform is preferable to overclaimed black-box scoring for early-stage prioritization

## Methods Sections

### Data Assembly

- collection and curation of viral protein-ligand structural records
- harmonization of virus, protein, ligand, and PDB-linked metadata

### Ligand Annotation Layers

- ligand synonym mapping
- atom-level ligand coordinate extraction
- SMILES association and atom-mapping workflows
- functional-group annotation

### Interaction and Surface Analysis

- protein-ligand contact derivation
- solvent-exposed atom assignment using SASA-derived workflows
- optional water or pocket-context layers where available

### Web Platform Implementation

- Flask-based application architecture
- SQLite-backed structural annotation store
- interactive query, comparison, and export surfaces

### PROTACability Heuristic Layer

- ligand-centered warhead linkability
- target lysine accessibility
- structural-priority aggregation
- geometry-informed degrader-readiness cues

## Results / Feature Demonstration Sections

- target-centric retrieval of viral protein-ligand structures
- ligand-first cross-structure mapping and interaction review
- solvent-exposed atom and functional-group interpretation
- multi-structure ligand comparison
- transparent PROTACability triage case studies
- downstream handoff into companion design tools

## Suggested Figures

1. Platform overview schematic showing structure-to-design workflow
2. Application map with major modules and evidence layers
3. Example ligand interaction and solvent-exposure view
4. Ligand comparison panel across multiple viral structures
5. PROTACability evidence-layer diagram with interpretation guardrails
6. Companion ecosystem figure linking V-LiSEMOD, Warhead Hunter, PROTAC Builder, and E3 Ligandalyzer

## Limitations

- dependent on available co-crystal structures and annotation quality
- heuristic PROTACability layers are not experimental degradation validation
- ligand mapping and solvent-exposure interpretation are limited by upstream data quality
- viral target coverage is constrained by source structure availability

## Future Work

- expanded viral target coverage
- additional structure-quality and confidence metadata
- stronger provenance tracking for generated annotations
- benchmarked case studies for warhead/linker hypothesis generation
- optional API formalization and lightweight test fixtures

## Suggested Keywords

viral structural bioinformatics, protein-ligand interactions, solvent-accessible surface area, warhead design, induced proximity, PROTACability, degrader readiness, antiviral discovery

## Citation Placeholder

Citation details to be added once manuscript, preprint, or software note metadata is finalized.
