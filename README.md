# V-LiSEMOD

## Viral Ligand Solvent-Exposed Moiety Database

Structure-guided viral protein-ligand analysis for solvent exposure review, interaction interpretation, ligand comparison, and transparent PROTACability-style triage.

[![Overview](https://img.shields.io/badge/Overview-Platform%20Summary-16a34a?style=for-the-badge)](#overview)
[![Capabilities](https://img.shields.io/badge/Capabilities-Core%20Modules-0891b2?style=for-the-badge)](#core-platform-capabilities)
[![Documentation](https://img.shields.io/badge/Documentation-docs%2F-ea580c?style=for-the-badge&logo=readthedocs)](#documentation)
[![PROTACability](https://img.shields.io/badge/PROTACability-Interpretation-7c3aed?style=for-the-badge)](#protacability-interpretation-note)
[![Deployment](https://img.shields.io/badge/Deployment-Quick%20Start-2563eb?style=for-the-badge)](#quick-start)
[![Companion Tools](https://img.shields.io/badge/Companion%20Tools-Ecosystem-0f766e?style=for-the-badge)](#companion-tool-ecosystem)

## Overview

V-LiSEMOD is a Flask-based structural bioinformatics platform for exploring curated viral protein-ligand co-crystal structures and the ligand-centered evidence that matters for follow-on design. It brings together viral protein metadata, ligand identity layers, atom-level interaction context, solvent-exposed atom analysis, functional-group annotations, and degrader-readiness heuristics in a single web workflow.

The public repository is now no-login by default. The web app does not require `users.db`, user accounts, or Flask-Login to access the scientific workflows.

The platform is intended for manuscript-ready structural interpretation, ligand prioritization, and hypothesis generation rather than black-box prediction.

## Scientific Motivation

Viral protein-ligand structures are rich sources of medicinal chemistry insight, but the design-relevant details are usually split across structure files, annotations, interaction outputs, solvent accessibility calculations, and custom downstream scripts. V-LiSEMOD was built to make those layers easier to inspect together so researchers can ask practical questions such as:

- Which viral targets have structure-resolved ligand contexts worth revisiting?
- Which ligand atoms appear solvent-exposed and potentially modifiable?
- Which contacts look important enough to preserve during linker growth?
- Which structures look more promising for degrader-oriented follow-up?

## Core Platform Capabilities

| Module | What it supports |
|---|---|
| Structure Explorer | Virus → PDB → ligand selection, ligand imagery, and PyMOL-oriented structure review |
| Protein Query | Target-centric filtering and export-ready dataset assembly |
| Ligand Indexer | Ligand-first lookup across mapped structural contexts |
| Ligand Comparison | Multi-structure interaction and atom-burden comparison |
| Solvent Exposure Analysis | Surface-exposed ligand atom review using SASA-derived data |
| Interaction Visualization | Contact summaries and ligand interaction chart generation |
| PROTACability Assessment | Transparent structural-priority and design-readiness triage |
| PROTAC Builder Handoff | External continuation into degrader design workflows |

## Conceptual Workflow

1. Start from a viral protein-ligand co-crystal structure.
2. Inspect target, ligand, chain, and residue context.
3. Review interactions, ligand atom exposure, and functional-group annotations.
4. Identify atoms that may tolerate modification while preserving key contacts.
5. Interpret PROTACability-style evidence as a prioritization aid.
6. Hand promising warhead contexts into downstream design tools when appropriate.

## Application Map

| Area | Role |
|---|---|
| Home / Structure Explorer | Main entry point for structure-specific exploration and PyMOL session generation |
| About | Scientific framing and project context |
| Protein Query | Virus/protein/ligand filtering and export workflows |
| Ligand Indexer | Ligand-first mapped-context lookup |
| Ligand Comparison | Cross-structure comparison of interaction behavior |
| PROTACability Assessment | Target, protein, structure, and chain-level triage views |
| Drug GPT / BioGPT | Optional local assistant module when enabled |

## PROTACability Interpretation Note

In V-LiSEMOD, PROTACability refers to transparent structural-priority and degrader-readiness heuristics assembled from ligand-centered and target-centered evidence layers. It does not mean experimentally validated degradation, productive ternary complex formation, or guaranteed medicinal chemistry tractability.

Preferred interpretation:

- Warhead linkability asks whether a bound ligand appears to expose plausible linker-attachment atoms.
- Target lysine accessibility asks whether exposed lysines exist as target-side cues.
- Structural priority and readiness layers support triage and hypothesis generation.

## Companion Tool Ecosystem

V-LiSEMOD is positioned as part of a broader induced-proximity and design workflow:

- [Warhead Hunter](https://warheadhunter.com): solvent-exposed atom and warhead-focused follow-up
- [PROTAC Builder](https://protacbuilder.com): downstream linker/recruiter/warhead assembly workflows
- [E3 Ligandalyzer](https://e3ligandalyzer.com): E3 ligase recruiter and ligase-context exploration

## Documentation

[![App Guide](https://img.shields.io/badge/docs-App%20Guide-2563eb?style=flat-square)](docs/APP_GUIDE.md)
[![Database](https://img.shields.io/badge/docs-Database-0f766e?style=flat-square)](docs/DATABASE.md)
[![PROTACability](https://img.shields.io/badge/docs-PROTACability-7c3aed?style=flat-square)](docs/PROTACABILITY.md)
[![Deployment](https://img.shields.io/badge/docs-Deployment-ea580c?style=flat-square)](docs/DEPLOYMENT.md)
[![Maintenance](https://img.shields.io/badge/docs-Maintenance-475569?style=flat-square)](docs/MAINTENANCE.md)
[![Developer Notes](https://img.shields.io/badge/docs-Developer%20Notes-1d4ed8?style=flat-square)](docs/DEVELOPER_NOTES.md)
[![Manuscript Outline](https://img.shields.io/badge/docs-Manuscript%20Outline-166534?style=flat-square)](docs/MANUSCRIPT_OUTLINE.md)

| Topic | Link |
|---|---|
| Documentation index | [docs/README.md](docs/README.md) |
| Application guide | [docs/APP_GUIDE.md](docs/APP_GUIDE.md) |
| Database layers | [docs/DATABASE.md](docs/DATABASE.md) |
| PROTACability interpretation | [docs/PROTACABILITY.md](docs/PROTACABILITY.md) |
| Deployment | [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md) |
| Maintenance | [docs/MAINTENANCE.md](docs/MAINTENANCE.md) |
| Developer notes | [docs/DEVELOPER_NOTES.md](docs/DEVELOPER_NOTES.md) |
| Manuscript framing | [docs/MANUSCRIPT_OUTLINE.md](docs/MANUSCRIPT_OUTLINE.md) |

## Repository Scope

Included in this repository:

- Flask application code and templates
- public documentation
- helper/import scripts
- lightweight project configuration and environment files

Not intended for public GitHub inclusion:

- local databases and user records
- private credentials or tokens
- large regenerated structure bundles and caches
- local-only model weights and checkpoints
- generated exports and other ephemeral outputs

Authentication note:

- `users.db` is no longer used by the application and should not be created, provisioned, or committed.
- `FLASK_SECRET_KEY` is still recommended because Flask may use anonymous browser sessions for temporary UI state.

## Quick Start

For practical setup, environment variables, and deployment notes, use [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md).

At a high level:

1. Install the Python dependencies listed in `requirements.txt` or the Conda environment in `environment.yml`.
2. Provide the required local data files and environment variables.
3. Run the app with `python app.py` for the current local default (`127.0.0.1:5003`) or use your chosen production entrypoint.

Remote data note:

- The simple V-LiSEMOD lookup routes can run against a separate RANDY API by setting `VLISMOD_DATA_BACKEND=randy` or `auto`.
- For production-style RANDY routing, prefer `VLISMOD_BACKUP_URL=https://randy.rove-vernier.ts.net/backup/vlismod`.
- RANDY should be configured with `VLISMOD_DB_PATH` and either `VLISMOD_API_TOKEN` or a shared RANDY backup token env var.
- V-LiSEMOD should be configured with `VLISMOD_BACKUP_URL`, `RANDY_API_TOKEN`, and `VLISMOD_DATA_BACKEND`.
- `RANDY_API_BASE_URL` is still supported for local/dev compatibility and will default to `/api/vlismod` if you provide only a host URL.
- If no RANDY configuration is provided, V-LiSEMOD continues to use local `viral_data.db` behavior.

## Project Status

V-LiSEMOD is an active research-oriented web platform and documentation cleanup is ongoing. The repository currently reflects a mix of production-facing Flask routes, data-enrichment scripts, optional local-LLM hooks, and manuscript-oriented scientific framing.

## Citation

Citation details to be added after manuscript or software record finalization.

## Contact

Please use the repository issue tracker or the project/lab contact route that accompanies the manuscript or deployment.

## Repository Description

V-LiSEMOD is a structural bioinformatics web platform for curated viral protein-ligand exploration, solvent-exposed moiety analysis, interaction review, and transparent degrader-readiness triage.

## Practical Takeaway

V-LiSEMOD helps researchers move from viral co-crystal structures to interpretable ligand-modification hypotheses without collapsing solvent exposure, contact preservation, and degrader-readiness questions into a single overclaimed prediction.
