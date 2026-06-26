# Manuscript Outline

## Working Title

V-LiSEMOD: a structure-guided viral protein-ligand exploration platform for solvent-exposed moiety analysis and transparent degrader-readiness triage

## Abstract Scaffold

Viral protein-ligand co-crystal structures contain useful medicinal chemistry and chemical biology signals, but those signals are often fragmented across structural files, ligand annotations, interaction outputs, solvent accessibility calculations, and downstream scripts. V-LiSEMOD is a Flask-based structural bioinformatics platform that brings curated viral protein metadata, ligand identity layers, atom-level interaction context, solvent-exposed atom analysis, functional-group annotations, ligand comparison workflows, and transparent PROTACability-style heuristics into a single interface. The platform supports target-centric query, ligand-first indexing, cross-structure comparison, PyMOL-oriented export, and degrader-readiness triage while keeping heuristic evidence separate from experimental degradation claims. V-LiSEMOD is best positioned as a hypothesis-generation and design-review resource for antiviral structural biology, medicinal chemistry, and induced-proximity-inspired workflows.

## Introduction Rationale

- viral structural datasets continue to expand, but design-relevant interpretation remains fragmented
- ligand optimization requires structure, contact, and solvent-exposure context in the same workflow
- degrader-oriented thinking adds target-side questions such as lysine accessibility and linker-attachment opportunity
- transparent evidence presentation is more defensible than overclaimed black-box prediction in early-stage triage

## Methods

### Software Implementation

- Flask-based web application
- route-heavy application architecture in `app.py`
- optional assistant blueprint in `DRUGapp.py`
- local SQLite mode and optional RANDY-backed route groups

### Database and Data Layers

- curated viral protein metadata
- ligand synonym and identifier mapping
- atom-level ligand coordinates
- contact layers and interaction summaries
- SASA-derived solvent-exposed atom data
- functional-group annotations
- PROTACability-linked tables imported from generated CSV layers

### App Modules

- Structure Explorer
- Protein Query
- Ligand Indexer
- Ligand Comparison
- PROTACability Assessment
- optional Drug GPT assistant surface

### PROTACability Heuristic Layer

- ligand-centered warhead linkability
- target-side lysine accessibility cues
- structural-priority aggregation
- combined degrader-readiness heuristic outputs

## Results or Feature Demonstration Structure

- target-centric retrieval of viral protein-ligand structures
- ligand-first mapped-context lookup and interaction review
- solvent-exposed atom and functional-group interpretation
- cross-structure ligand comparison
- transparent PROTACability triage example workflows
- companion-tool handoff into downstream design environments

## Limitations

- dependent on available co-crystal structures and annotation quality
- dependent on local or provisioned data availability
- PROTACability is heuristic and not experimental degradation validation
- geometry and lysine-context cues do not guarantee productive ternary-complex formation
- optional assistant behavior is deployment-dependent
- generated outputs depend on writable local runtime storage

## Future Work

- broader viral target coverage
- stronger provenance tracking and fixture data for reproducibility
- benchmarking and structured case studies
- more formal API and test coverage
- expanded companion-tool interoperability

## Figure Ideas

1. Connected ecosystem overview including V-LiSEMOD, Warhead Hunter, PROTAC Builder, E3 Ligandalyzer, and PyMACS
2. Application map covering the main V-LiSEMOD user workflows
3. Warhead Hunter workflow and results-library figure
4. V-LiSEMOD workflow composite with Structure Explorer, Protein Query, Ligand Indexer, Ligand Comparison, and PROTACability
5. PROTAC Builder downstream continuation figure
6. E3 Ligandalyzer and PyMACS companion-context figure
7. Public API or reproducibility evidence figure using lightweight public endpoint checks

## Keywords

viral structural bioinformatics, protein-ligand interactions, solvent-accessible surface area, ligand modification, induced proximity, PROTACability, degrader readiness, antiviral discovery

## Manuscript-Safe Language

- “V-LiSEMOD is a Flask-based structural bioinformatics platform for exploring curated viral protein-ligand co-crystal contexts.”
- “The platform integrates ligand identity, atom-level interaction context, solvent exposure evidence, functional-group annotations, and transparent degrader-readiness heuristics.”
- “PROTACability outputs are intended for hypothesis generation and design triage rather than experimental degradation prediction.”
- “Companion-tool links support continuity with downstream design workflows, including PROTAC Builder and related structural-analysis tools.”

## Language To Avoid Unless Later Validated

- “validated degrader prediction”
- “guaranteed PROTACability”
- “automated PROTAC design”
- “complete medicinal chemistry decision engine”
- “exhaustive viral ligand database”
- “production-hardened public API”
- “experimentally confirmed degradation”

## Citation Placeholder

Add manuscript, preprint, or software-record citation details after submission or release.
