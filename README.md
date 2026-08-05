# V-LiSEMOD

<p align="center">
  <strong>Viral Ligand Structure Explorer for Modification-Oriented Design</strong>
</p>

<p align="center">
  <em>A Flask-based structural bioinformatics platform for exploring viral protein-ligand co-crystals, ligand interaction evidence, solvent-exposed atoms, and transparent degrader-readiness triage.</em>
</p>

<p align="center">
  <a href="https://vlisemod.com">
    <img src="https://img.shields.io/badge/Open-V--LiSEMOD-1f6feb?style=for-the-badge" alt="Open V-LiSEMOD">
  </a>
  <a href="https://vlisemod.com/protacability_page">
    <img src="https://img.shields.io/badge/Open-PROTACability-2ea44f?style=for-the-badge" alt="Open PROTACability assessment">
  </a>
  <a href="https://vlisemod.com/query_protein_virus_page">
    <img src="https://img.shields.io/badge/Open-Protein%20Query-f97316?style=for-the-badge" alt="Open Protein Query">
  </a>
  <a href="https://vlisemod.com/contact">
    <img src="https://img.shields.io/badge/Contact-Collaboration-6f42c1?style=for-the-badge" alt="Contact V-LiSEMOD">
  </a>
</p>

---

## Overview

V-LiSEMOD helps researchers move from curated viral protein-ligand structures to interpretable ligand-modification hypotheses. It brings together structure selection, ligand interaction context, solvent exposure review, functional-group annotations, and PROTACability-style triage in one web application.

The project is intended for structural review, hypothesis generation, manuscript support, and collaborative tool development. It is not a validated degradation-prediction engine and should not be interpreted as a guarantee of PROTAC design success.

## Live Tool

The public V-LiSEMOD entry point is:

**[https://vlisemod.com](https://vlisemod.com)**

Direct pages:

| Page | Link | Purpose |
|---|---|---|
| Structure Explorer | [vlisemod.com](https://vlisemod.com) | Select virus, PDB structure, ligand, chain, and residue context for structure-centered review. |
| About | [vlisemod.com/about](https://vlisemod.com/about) | Project framing, scientific motivation, workflow overview, and ecosystem context. |
| Protein Query | [vlisemod.com/query_protein_virus_page](https://vlisemod.com/query_protein_virus_page) | Filter by virus, target, structure, and ligand; assemble export-oriented datasets. |
| Ligand Indexer | [vlisemod.com/ligand_indexer](https://vlisemod.com/ligand_indexer) | Search ligands and synonyms across mapped viral co-crystal contexts. |
| Ligand Comparison | [vlisemod.com/compare_ligands](https://vlisemod.com/compare_ligands) | Compare ligand interaction behavior across multiple structural contexts. |
| PROTACability Assessment | [vlisemod.com/protacability_page](https://vlisemod.com/protacability_page) | Review transparent structural-priority and degrader-readiness heuristics. |
| Use Cases | [vlisemod.com/use-cases](https://vlisemod.com/use-cases) | Follow common research workflows through the application. |
| Viral PROTAC Design | [vlisemod.com/viral-protac-design](https://vlisemod.com/viral-protac-design) | Read the design rationale behind viral degrader-oriented triage. |
| Viral Drug Targets | [vlisemod.com/viral-drug-targets](https://vlisemod.com/viral-drug-targets) | Review broader target-context material for structure-guided antiviral design. |
| In Silico Virology Tools | [vlisemod.com/in-silico-virology-tools](https://vlisemod.com/in-silico-virology-tools) | Explore the companion-tool ecosystem. |
| Methods | [vlisemod.com/methods](https://vlisemod.com/methods) | Review method notes and interpretation guidance. |
| FAQ | [vlisemod.com/faq](https://vlisemod.com/faq) | Answer common interpretation and usage questions. |
| Citation | [vlisemod.com/citation](https://vlisemod.com/citation) | Find citation and attribution guidance. |
| Contact | [vlisemod.com/contact](https://vlisemod.com/contact) | Start a collaboration or methods discussion. |

## Core Capabilities

| Module | What it supports |
|---|---|
| Structure Explorer | Virus -> PDB -> ligand exploration, ligand imagery, PyMOL-oriented exports, solvent exposure review, and functional-group context. |
| Protein Query | Target-centric filtering, PDB selection, ligand-aware lookup, and export-oriented dataset assembly. |
| Ligand Indexer | Ligand-first search across mapped viral protein-ligand structure contexts. |
| Ligand Comparison | Multi-structure comparison of ligand interaction behavior and atom-level contact patterns. |
| Solvent Exposure Review | SASA-derived exposed-atom interpretation for modification-site discussion. |
| PROTACability Assessment | Transparent target, structure, chain, lysine-proximity, warhead-linkability, and readiness triage. |
| Companion Tool Handoff | Continuation into external design environments such as PROTAC Builder. |

## Conceptual Workflow

1. Start from a viral protein-ligand co-crystal structure.
2. Inspect virus, target, PDB, ligand, chain, and residue context.
3. Review interaction evidence, solvent-exposed atoms, and functional-group annotations.
4. Identify ligand-centered modification opportunities while preserving important contacts.
5. Interpret PROTACability results as triage cues for follow-up design review.
6. Continue promising contexts in downstream companion tools when appropriate.

## PROTACability Interpretation

In V-LiSEMOD, PROTACability means transparent structural-priority and design-readiness heuristics assembled from ligand-centered and target-centered evidence layers.

Use these outputs as triage signals:

- warhead linkability highlights plausible ligand-centered attachment opportunities,
- lysine accessibility provides target-side structural cues,
- structure priority helps order follow-up review,
- degrader-readiness layers summarize heuristic context.

Do not describe the output as guaranteed PROTAC design success, productive ternary-complex prediction, or experimentally validated viral protein degradability.

## Companion Tool Ecosystem

V-LiSEMOD sits upstream of a broader structure-guided induced-proximity workflow:

| Tool | Link | Role |
|---|---|---|
| Warhead Hunter | [warheadhunter.com](https://warheadhunter.com) | Warhead-focused follow-up and solvent-exposed atom context. |
| PROTAC Builder | [protacbuilder.com](https://protacbuilder.com) | Downstream linker, recruiter, and bifunctional design workflows. |
| E3 Ligandalyzer | [e3ligandalyzer.com](https://e3ligandalyzer.com) | E3 ligase and recruiter-context exploration. |

These tools are conceptually connected. They are not all implemented inside this repository.

## Repository Scope

This repository is prepared for public source review. It includes:

- Flask application code,
- templates and static assets,
- documentation for application behavior, interpretation, development, and maintenance,
- enrichment, import, and maintenance scripts,
- small curated source assets needed to understand the project.

It intentionally excludes private runtime assets such as:

- credentials, tokens, and private configuration,
- generated structure, chart, export, and cache directories,
- generated model, training, and runtime artifacts,
- local backups and archives.

## Development Note

The public repository is shared for transparency, review, and collaboration around the V-LiSEMOD application code and scientific framing. The hosted application at [vlisemod.com](https://vlisemod.com) is the primary way to use the tool.

Local development is possible for collaborators with project access, but this README focuses on the public tool rather than setup instructions.

## Documentation

- [Documentation index](docs/README.md)
- [Application guide](docs/APP_GUIDE.md)
- [PROTACability guide](docs/PROTACABILITY.md)
- [Developer notes](docs/DEVELOPER_NOTES.md)
- [Maintenance guide](docs/MAINTENANCE.md)

## Public Release Checklist

Before moving this repository from private to public:

- confirm local credential and configuration files are untracked,
- confirm private runtime data, backups, generated caches, and model artifacts are untracked,
- run a secret scan over tracked files,
- verify documentation links point to tracked files or live public pages,
- verify public claims remain aligned with the implementation and live tool behavior.

Helpful commands:

```bash
git status --ignored --short
git ls-files -- '*secret*' '*credentials*' '*.pem' '*.key'
rg -n "(api[_-]?key|secret|token|password|credential|private[_-]?key|sk-[A-Za-z0-9]|ghp_|github_pat_|AKIA[0-9A-Z]{16})" --glob "!.git/**"
```

## Project Status

V-LiSEMOD is an active research-oriented web platform. The public repository supports review, collaboration, and development around the application code, documentation, and scientific framing. The hosted site remains the primary public tool.

## Citation

Citation details can be added once a manuscript, preprint, or software record is finalized. Until then, cite the V-LiSEMOD website and repository according to the guidance on [vlisemod.com/citation](https://vlisemod.com/citation).

## Practical Takeaway

V-LiSEMOD helps users move from curated viral co-crystal structures to interpretable ligand-modification and degrader-follow-up hypotheses without hiding the decision process inside a black-box predictor.
