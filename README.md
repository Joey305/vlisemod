# V-LiSEMOD

V-LiSEMOD is a Flask-based structural bioinformatics platform for exploring curated viral protein-ligand co-crystal contexts, reviewing ligand-centered interaction evidence, and triaging solvent-exposed modification opportunities with transparent PROTACability-style heuristics.

The repository is intended to support reviewers, collaborators, manuscript preparation, and future development. It describes a research-oriented web application, not a validated degradation-prediction engine.

## Scientific Motivation

Viral protein-ligand structures contain practical design signals, but those signals are often split across structure files, atom-level contact tables, solvent accessibility workflows, and downstream scripts. V-LiSEMOD brings those evidence layers into one interface so users can ask:

- which viral targets have useful ligand-bound structural context,
- which ligand atoms appear solvent-exposed and potentially modifiable,
- which contacts may need to be preserved during optimization, and
- which structure contexts may merit degrader-oriented follow-up.

## Core Capabilities

| Module | Current role |
|---|---|
| Structure Explorer | Virus -> PDB -> ligand exploration, ligand imagery, and PyMOL-oriented structural review |
| Protein Query | Target-centric filtering, structure selection, and export-oriented dataset assembly |
| Ligand Indexer | Ligand-first lookup across mapped viral structure contexts |
| Ligand Comparison | Cross-structure comparison of ligand interaction behavior and mapped context |
| Solvent Exposure Review | SASA-derived exposed-atom interpretation and functional-group context |
| PROTACability Assessment | Transparent structural-priority and degrader-readiness triage |
| Companion-tool handoff | External continuation into PROTAC Builder and related tools |
| Optional Drug GPT module | Deployment-dependent assistant workflow when feature flags and local model runtime are enabled |

## Conceptual Workflow

1. Start from a viral protein-ligand co-crystal structure.
2. Inspect target, ligand, chain, and residue context.
3. Review interaction evidence, solvent-exposed atoms, and functional-group annotations.
4. Identify ligand-centered modification opportunities while preserving important contacts.
5. Interpret PROTACability outputs as triage cues for follow-up design review.
6. Hand promising contexts into downstream companion tools when appropriate.

## Application Map

| Area | What it supports |
|---|---|
| Home / Structure Explorer | Structure-specific exploration and export-oriented review |
| About | Project framing, workflow explanation, and ecosystem positioning |
| Protein Query | Virus/protein/ligand filtering and data-export assembly |
| Ligand Indexer | Ligand-centered mapped-context lookup |
| Ligand Comparison | Multi-structure interaction and context comparison |
| PROTACability Assessment | Target, protein, structure, and chain-level heuristic triage |
| Drug GPT / BioGPT | Optional assistant surface when explicitly enabled |

## PROTACability Interpretation Note

In V-LiSEMOD, PROTACability refers to transparent structural-priority and design-readiness heuristics assembled from ligand-centered and target-centered evidence layers. It is intended for hypothesis generation and triage rather than experimental degradation prediction.

Safe present-tense interpretation:

- warhead linkability highlights plausible ligand-centered attachment opportunity,
- lysine accessibility provides target-side structural cues,
- structural priority supports follow-up review ordering, and
- combined readiness layers summarize heuristic degrader-oriented context.

It should not be described as guaranteed PROTAC design success, productive ternary-complex prediction, or experimentally validated degradability.

## Companion Tool Ecosystem

V-LiSEMOD is framed as part of a broader structure-guided induced-proximity workflow:

- Warhead Hunter for warhead-focused follow-up and solvent-exposed atom context
- PROTAC Builder for downstream linker and recruiter design workflows
- E3 Ligandalyzer for E3 ligase and recruiter-context exploration

These companion tools are conceptually connected workflows. They are not fully implemented inside this repository.

## Repository Scope

Included here:

- Flask application code
- templates and static assets
- documentation and manuscript-planning materials
- helper scripts for enrichment, import, and maintenance

Not intended for public GitHub inclusion:

- local databases such as `viral_data.db`
- generated caches, exports, and PyMOL sessions
- local credentials or tokens
- model weights or checkpoints
- large regenerated archives and runtime downloads

## Quick Start

1. Install dependencies from `requirements.txt` or `environment.yml`.
2. Provide a local `viral_data.db` or configure RANDY-backed data access.
3. Set the environment variables needed for your deployment mode.
4. Start the app with `python app.py` for local development.
5. Validate `/healthz` and the main public pages.

The app supports both a local SQLite-backed mode and a RANDY-backed mode for selected route groups. The optional Drug GPT module remains deployment-dependent and is disabled by default unless enabled explicitly.

## Documentation

- [Documentation index](docs/README.md)
- [Application guide](docs/APP_GUIDE.md)
- [Database guide](docs/DATABASE.md)
- [PROTACability guide](docs/PROTACABILITY.md)
- [Deployment guide](docs/DEPLOYMENT.md)
- [Developer notes](docs/DEVELOPER_NOTES.md)
- [Maintenance guide](docs/MAINTENANCE.md)
- [Manuscript outline](docs/MANUSCRIPT_OUTLINE.md)
- [Detailed manuscript planning](docs/manuscript/MANUSCRIPT_PLAN.md)

## Project Status

V-LiSEMOD is an active research-oriented web platform. The repository currently combines a live Flask application, local or provisioned data dependencies, optional companion-tool integrations, and manuscript-planning material intended to keep public claims aligned with the actual implementation.

## Citation

Citation details can be added once a manuscript, preprint, or software record is finalized.

## Practical Takeaway

V-LiSEMOD helps users move from curated viral co-crystal structures to interpretable ligand-modification and degrader-follow-up hypotheses without collapsing those decisions into an overclaimed black-box predictor.
