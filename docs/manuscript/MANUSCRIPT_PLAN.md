# Manuscript Plan

## Working Title Options

1. V-LiSEMOD: a web platform for viral protein-ligand exploration and solvent-exposed moiety analysis
2. V-LiSEMOD: structure-guided viral ligand analysis with transparent degrader-readiness triage
3. V-LiSEMOD: a Flask-based viral structural bioinformatics resource for ligand-centered design review
4. V-LiSEMOD: a curated viral protein-ligand exploration platform for solvent exposure and PROTACability-style assessment
5. V-LiSEMOD: integrating viral ligand context, solvent exposure, and degrader-oriented heuristic triage
6. V-LiSEMOD: a research web application for viral co-crystal interpretation and ligand comparison
7. V-LiSEMOD: a structure-guided resource for viral ligand interrogation and design-readiness review

## Manuscript Thesis

V-LiSEMOD is best framed as a structure-guided viral protein-ligand exploration platform that integrates curated target metadata, ligand-centered interaction evidence, solvent-exposed atom review, functional-group annotations, cross-structure comparison, and transparent degrader-readiness heuristics into a single web workflow. Its contribution is not a claim of validated degradation prediction, but a practical interface for hypothesis generation, design triage, and companion-tool handoff in antiviral and induced-proximity-inspired research.

## Target Article Type

- software article
- webserver paper
- application note
- methods or resource paper

Do not claim fit for a specific journal without checking that journal’s current scope and submission requirements.

## Intended Audience

- structural biologists
- medicinal chemists
- chemical biologists
- antiviral researchers
- PROTAC and degrader researchers
- computational chemists
- developers extending structural-analysis tools

## Core Contributions

- curated viral protein-ligand exploration
- ligand-centered interaction review
- solvent-exposed atom interpretation
- functional-group and atom-level context
- ligand comparison across structures
- transparent PROTACability-style heuristic triage
- companion-tool handoff into downstream design workflows

## Manuscript Positioning

- V-LiSEMOD as a web platform for viral structure-linked ligand interpretation
- V-LiSEMOD as a database-backed analysis interface
- V-LiSEMOD as a hypothesis-generation and triage tool
- V-LiSEMOD as one component of a broader companion-tool ecosystem

## Claims To Emphasize

- Flask-based web application
- public no-login workflow by default in the current app
- SQLite-backed curated data layers in local mode
- optional RANDY-backed route groups when configured
- Structure Explorer workflow
- Protein Query workflow
- Ligand Indexer workflow
- Ligand Comparison workflow
- PROTACability Assessment workflow
- companion-tool handoff to PROTAC Builder
- optional local assistant module only when enabled

## Claims To Avoid

- experimentally validated degradation prediction
- guaranteed PROTAC design success
- productive ternary-complex prediction
- fully automated medicinal chemistry decision-making
- exhaustive viral ligand coverage
- production-hardened hosted database infrastructure
- authentication or rate-limiting claims unless later implemented and documented
- full formal API claims unless later documented and validated
