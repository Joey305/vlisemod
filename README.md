# 🧬 V-LiSEMOD

<p align="center">
  <strong>Viral Ligand Solvent-Exposed Moiety Database</strong>
</p>

<p align="center">
  <em>A structure-guided web platform for exploring viral protein–ligand co-crystal structures, ligand interaction networks, solvent-exposed ligand atoms, warhead/linkability evidence, and degrader-readiness triage.</em>
</p>

<p align="center">
  <a href="#overview">
    <img src="https://img.shields.io/badge/Project-Overview-22c55e?style=for-the-badge" alt="Project Overview">
  </a>
  <a href="#platform-capabilities">
    <img src="https://img.shields.io/badge/Explore-Capabilities-06b6d4?style=for-the-badge" alt="Capabilities">
  </a>
  <a href="#scientific-framing">
    <img src="https://img.shields.io/badge/Scientific-Framing-blueviolet?style=for-the-badge" alt="Scientific Framing">
  </a>
  <a href="#documentation">
    <img src="https://img.shields.io/badge/Docs-Project%20Documentation-orange?style=for-the-badge&logo=readthedocs" alt="Documentation">
  </a>
</p>

<p align="center">
  <a href="https://warheadhunter.com">
    <img src="https://img.shields.io/badge/Companion-Warhead%20Hunter-00e5ff?style=for-the-badge" alt="Warhead Hunter">
  </a>
  <a href="https://protacbuilder.com">
    <img src="https://img.shields.io/badge/Companion-PROTAC%20Builder-06b6d4?style=for-the-badge" alt="PROTAC Builder">
  </a>
  <a href="https://e3ligandalyzer.com">
    <img src="https://img.shields.io/badge/Companion-E3%20Ligandalyzer-7c3aed?style=for-the-badge" alt="E3 Ligandalyzer">
  </a>
</p>

---

<p align="center">
  <strong>From viral protein–ligand structures to interpretable warhead and degrader-design evidence.</strong>
</p>

<p align="center">
  <em>V-LiSEMOD connects structural biology, cheminformatics, ligand interaction analysis, and targeted protein degradation workflows.</em>
</p>

---

<a id="overview"></a>

## 🚀 Overview

**V-LiSEMOD** — the **Viral Ligand Solvent-Exposed Moiety Database** — is a structural bioinformatics platform for exploring curated viral protein–ligand co-crystal structures and extracting design-relevant ligand evidence.

The platform was built to help researchers answer questions such as:

```text
Which viral protein-ligand structures are available for a target?
How does a ligand interact across different viral structures?
Which ligand atoms are solvent-exposed and potentially modifiable?
Which atoms appear interaction-critical and should be protected?
Which structures may be useful starting points for degrader-oriented design?
```

V-LiSEMOD combines viral target metadata, ligand identifiers, protein–ligand interaction data, solvent accessibility measurements, ligand atom annotations, SMILES-linked mappings, and PROTACability-style evidence layers into a searchable web interface.

The project is designed as part of a broader induced-proximity and molecular design ecosystem that includes **Warhead Hunter**, **PROTAC Builder**, and **E3 Ligandalyzer**.

---

<a id="scientific-motivation"></a>

## 🎯 Scientific Motivation

Viral protein–ligand co-crystal structures contain valuable information for antiviral inhibitor discovery, covalent warhead design, linker placement, and targeted degrader development.

However, the information needed for structure-guided design is rarely available in one place. Researchers often need to manually connect:

```text
PDB structures
ligand names and synonyms
protein target annotations
ligand atom identities
protein-ligand contact networks
solvent-accessible surface area
functional group annotations
SMILES / atom mappings
potential linker attachment atoms
target lysine accessibility
```

V-LiSEMOD brings these layers together so that ligand modification hypotheses can be reviewed in the context of the actual bound structure.

---

<a id="platform-capabilities"></a>

## ✨ Platform Capabilities

| Module | Purpose |
|---|---|
| **Structure Explorer** | Search and inspect viral protein–ligand structures by virus, PDB, ligand, and chain context. |
| **Protein Query** | Filter viral structures by target, protein class, ligand, and available annotation layers. |
| **Ligand Indexer** | Start from a ligand code or synonym and find mapped structural contexts. |
| **Ligand Comparison** | Compare one ligand across structures using interaction distributions, distance profiles, and atom-level burden. |
| **Solvent-Exposed Atom Analysis** | Identify ligand atoms that may tolerate chemical modification or linker growth. |
| **Interaction Diagramming** | Visualize ligand-protein interaction patterns and contact distributions. |
| **PyMOL / Structure Export** | Generate structure-viewing outputs for downstream manual inspection. |
| **PROTACability Assessment** | Triage viral target/structure/ligand combinations using transparent degrader-readiness evidence layers. |
| **PROTAC Builder Handoff** | Move promising ligand contexts into downstream degrader design workflows. |

---

<a id="conceptual-workflow"></a>

## 🧩 Conceptual Workflow

V-LiSEMOD is organized around a structure-to-design workflow:

```text
Viral protein-ligand co-crystal structure
        ↓
Target, ligand, chain, and residue context
        ↓
Protein-ligand interaction analysis
        ↓
Ligand atom solvent exposure
        ↓
Functional group and SMILES-linked annotations
        ↓
Candidate modifiable atom / warhead evidence
        ↓
Target lysine accessibility and structural-priority evidence
        ↓
PROTACability-style triage
        ↓
Downstream design handoff
```

The goal is not to replace expert review, but to make the structural evidence easier to find, compare, export, and interpret.

---

<a id="application-map"></a>

## 🧭 Application Map

| Page / Area | What it supports |
|---|---|
| **Home / Structure Explorer** | Structure-specific exploration, ligand views, PyMOL outputs, and solvent-exposed atom visualization. |
| **About** | Scientific context, project motivation, and ecosystem positioning. |
| **Protein Query** | Target-centric search and export workflows. |
| **Ligand Indexer** | Ligand-first lookup across PDB, chain, and residue contexts. |
| **Ligand Comparison** | Multi-structure ligand interaction and atom-burden comparison. |
| **PROTACability** | Degrader-readiness triage using separable evidence layers. |
| **PROTAC Builder Linkout** | Downstream design handoff for linker/recruiter exploration. |

Detailed route, template, database, and maintenance notes are intentionally kept outside the public README. See [Documentation](#documentation).

---

<a id="protacability"></a>

## 🧪 PROTACability Assessment

The PROTACability module organizes degrader-relevant evidence into interpretable layers rather than treating the result as a black-box prediction.

| Evidence layer | Interpretation |
|---|---|
| **Warhead Linkability** | Ligand-centered evidence that a bound ligand may contain solvent-exposed, chemically modifiable atoms suitable for linker attachment. |
| **Target Lysine Accessibility** | Target-side evidence that accessible lysines are present as ubiquitination-readiness cues. |
| **Protein Structural Priority** | Structure-level prioritization based on ligand context, protein context, lysine availability, pI-related cues, and structural evidence. |
| **Ternary Geometry Cue** | Hypothesis-generating cue based on the relationship between ligand context and exposed lysine proximity. |
| **Overall Degrader Readiness** | Combined triage score summarizing the available evidence layers. |

> **Important interpretation note:** PROTACability outputs are transparent structural-priority and design-readiness heuristics. They support hypothesis generation and triage, but they are not experimentally validated degradation predictions.

---

<a id="scientific-framing"></a>

## 🧠 Scientific Framing

V-LiSEMOD should be interpreted as a **structure-guided evidence platform**, not as an automatic medicinal chemistry decision-maker.

The platform is intended to support:

- antiviral structure-guided discovery,
- ligand modification analysis,
- warhead and exit-vector prioritization,
- cross-structure ligand comparison,
- degrader-readiness triage,
- manuscript-ready structural interpretation, and
- downstream design handoff.

Preferred language:

```text
V-LiSEMOD identifies structure-supported ligand modification opportunities.
V-LiSEMOD highlights solvent-exposed ligand atoms that may support linker attachment.
V-LiSEMOD separates ligand-centered warhead evidence from target-side lysine accessibility.
V-LiSEMOD provides hypothesis-generating degrader-readiness triage.
```

Avoid overclaiming:

```text
V-LiSEMOD does not experimentally validate degradation.
V-LiSEMOD does not prove productive ternary-complex formation.
V-LiSEMOD does not replace medicinal chemistry review or experimental validation.
```

---

<a id="data-scope"></a>

## 🧾 Data Scope

V-LiSEMOD is built around curated and derived structural data layers, including:

- viral protein-ligand co-crystal structures,
- ligand and synonym mappings,
- ligand atom records,
- protein-ligand interaction/contact data,
- solvent-accessibility measurements,
- functional group annotations,
- SMILES and PDB-to-SMILES atom mappings,
- receptor binding-pocket annotations,
- PROTACability assessment tables,
- warhead/linkability enrichment data,
- target lysine accessibility data.

Large structural files, generated outputs, local databases, and private/internal datasets are not expected to live directly in the public source repository unless a specific lightweight example bundle is provided.

---

<a id="companion-tool-ecosystem"></a>

## 🧬 Companion Tool Ecosystem

V-LiSEMOD is part of a broader modular workflow for structure-guided molecular design and induced-proximity research.

```text
V-LiSEMOD
    Viral protein-ligand structure exploration
    Ligand interactions
    Solvent-exposed moiety analysis
    PROTACability triage

Warhead Hunter
    Structure-aware solvent-exposed atom analysis
    Exit-vector and warhead modification support

PROTAC Builder
    Linker / recruiter / warhead assembly
    PROTAC-like molecule construction
    Descriptor and reporting workflows

E3 Ligandalyzer
    E3 recruiter and ligase-focused structure analytics
```

Together, these tools support a connected workflow:

```text
protein-ligand structure
        ↓
modifiable atom / exit-vector evidence
        ↓
E3 recruiter and linker exploration
        ↓
PROTAC-like molecule generation
        ↓
modeling, simulation, and experimental prioritization
```

---

<a id="documentation"></a>

## 📚 Documentation

This README is intentionally high-level.

Detailed setup, maintenance, deployment, and handoff instructions should live in dedicated documentation files, for example:

| Document | Purpose |
|---|---|
| `docs/APP_GUIDE.md` | Page-by-page user guide and feature walkthrough. |
| `docs/DATABASE.md` | Database schema, table descriptions, and data-layer notes. |
| `docs/PROTACABILITY.md` | PROTACability scoring logic, interpretation, and regeneration notes. |
| `docs/DEPLOYMENT.md` | Deployment steps for Heroku, Azure, or other Flask-compatible platforms. |
| `docs/MAINTENANCE.md` | Local maintenance, cache cleanup, generated output policy, and troubleshooting. |
| `docs/DEVELOPER_NOTES.md` | Route inventory, optional modules, environment variables, and implementation notes. |
| `docs/MANUSCRIPT_OUTLINE.md` | Manuscript framing, figure ideas, limitations, and future work. |

This keeps the public README readable while preserving the full technical handoff for collaborators and maintainers.

---

<a id="repository-scope"></a>

## 📦 Repository Scope

The public-facing repository is intended to contain:

- application source code,
- templates and static assets,
- documentation,
- lightweight example or demo assets where appropriate,
- reproducible configuration examples,
- route and feature descriptions,
- manuscript-oriented documentation.

The repository is not intended to contain:

- large PDB structure folders,
- local SQLite databases with full internal datasets,
- generated PyMOL sessions,
- generated ligand image caches,
- generated charts,
- local `.env` files,
- API keys or tokens,
- model checkpoints,
- private/internal lab data.

A complete deployment may require external data files or database assets described in the documentation.

---

<a id="quick-start"></a>

## ⚡ Quick Start

This repository is organized as a Flask application.

For installation, environment setup, database configuration, and deployment instructions, see:

```text
docs/DEPLOYMENT.md
docs/MAINTENANCE.md
docs/DATABASE.md
```

A minimal local development flow is expected to follow the standard pattern:

```bash
git clone <repository-url>
cd VLISEMOD
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```

The exact runtime configuration depends on the available database, generated assets, and deployment environment.

---

<a id="status"></a>

## 🚧 Project Status

V-LiSEMOD is under active development.

Current priorities include:

- improving public-facing documentation,
- separating public source code from private/generated data assets,
- polishing light/dark UI consistency,
- improving graph and comparison layouts,
- refining PROTACability interpretation and reporting,
- strengthening route/deployment validation,
- preparing manuscript-oriented figures and documentation,
- improving downstream PROTAC Builder integration.

---

<a id="citation"></a>

## 🧬 Citation

A formal manuscript/software citation will be added when available.

Suggested interim citation:

```text
Schulz, J.-M. V-LiSEMOD: Viral Ligand Solvent-Exposed Moiety Database for Structure-Guided Warhead Discovery and Degrader-Readiness Triage. Software platform and repository, in development.
```

Related work motivating the degrader-readiness framing includes:

```text
Khurshid et al. Targeted degrader technologies as prospective SARS-CoV-2 therapies. Drug Discovery Today. 2024.
```

---

<a id="contact"></a>

## 📬 Contact

For collaboration, questions, or project discussion:

<p align="center">
  <a href="mailto:jxs794@miami.edu?subject=V-LiSEMOD%20Question%20%2F%20Collaboration">
    <img src="https://img.shields.io/badge/Joseph--Michael%20Schulz-jxs794%40miami.edu-blue?style=for-the-badge&logo=gmail" alt="Email Joseph-Michael Schulz">
  </a>
</p>

---

<a id="repository-description"></a>

## 🧾 Repository Description

> V-LiSEMOD is a structure-guided viral ligand database and web platform for exploring protein-ligand interactions, solvent-exposed ligand atoms, warhead/linkability evidence, and PROTACability-style degrader-readiness triage.

---

<a id="practical-takeaway"></a>

## 🙌 Practical Takeaway

Use V-LiSEMOD to move from:

```text
viral protein-ligand structures
```

to:

```text
interpretable ligand interaction evidence,
solvent-exposed modification opportunities,
warhead/linkability hypotheses,
and degrader-readiness triage.
```

The platform is designed to make viral structural data more useful for rational inhibitor optimization, warhead discovery, and induced-proximity design.
