# 🧬 V-LiSEMOD

<p align="center">
  <strong>V-LiSEMOD: Viral Ligand Solvent-Exposed Moiety Database for Structure-Guided Warhead Discovery, Ligand Interaction Analysis, and PROTACability Triage</strong>
</p>

<p align="center">
  <em>A Flask-based structural bioinformatics platform for exploring viral protein–ligand co-crystal structures, solvent-exposed ligand atoms, protein–ligand interactions, warhead linkability evidence, and degrader-readiness heuristics.</em>
</p>

<p align="center">
  <a href="#quick-start">
    <img src="https://img.shields.io/badge/Get%20Started-Quick%20Start-orange?style=for-the-badge&logo=gnubash" alt="Quick Start">
  </a>
  <a href="#core-capabilities">
    <img src="https://img.shields.io/badge/Explore-Core%20Capabilities-22c55e?style=for-the-badge" alt="Core Capabilities">
  </a>
  <a href="#application-pages">
    <img src="https://img.shields.io/badge/App%20Guide-Page%20Map-06b6d4?style=for-the-badge" alt="Application Pages">
  </a>
  <a href="#deployment">
    <img src="https://img.shields.io/badge/Deployment-GitHub%20%2B%20Heroku-blueviolet?style=for-the-badge&logo=heroku" alt="Deployment">
  </a>
</p>

<p align="center">
  <a href="https://protacbuilder.com">
    <img src="https://img.shields.io/badge/Companion%20Tool-PROTAC%20Builder-06b6d4?style=for-the-badge" alt="PROTAC Builder">
  </a>
  <a href="https://warheadhunter.com">
    <img src="https://img.shields.io/badge/Companion%20Tool-Warhead%20Hunter-00e5ff?style=for-the-badge" alt="Warhead Hunter">
  </a>
  <a href="https://e3ligandalyzer.com">
    <img src="https://img.shields.io/badge/Companion%20Tool-E3%20Ligandalyzer-7c3aed?style=for-the-badge" alt="E3 Ligandalyzer">
  </a>
</p>

<p align="center">
  <a href="mailto:jxs794@miami.edu?subject=V-LiSEMOD%20Question%20%2F%20Collaboration">
    <img src="https://img.shields.io/badge/Contact-Joseph--Michael%20Schulz-blue?style=for-the-badge&logo=gmail" alt="Contact Joseph-Michael Schulz">
  </a>
</p>

---

<p align="center">
  <strong>Search viral structures. Inspect ligand contacts. Find solvent-exposed atoms. Prioritize warhead/linker opportunities.</strong>
</p>

<p align="center">
  <em>From viral protein–ligand co-crystal data to interpretable ligand modification evidence and downstream degrader design handoff.</em>
</p>

---

<a id="overview"></a>

## 🚀 Overview

**V-LiSEMOD** — the **Viral Ligand Solvent-Exposed Moiety Database** — is a web-based platform for exploring curated viral protein–ligand structures from the Protein Data Bank and related structural/chemical annotations.

The application helps researchers move from broad target or ligand questions to detailed structural evidence, including:

- viral protein and ligand search,
- protein-centric query workflows,
- ligand-first indexing,
- protein–ligand interaction analysis,
- atom-level ligand contact burden,
- solvent-exposed ligand atom visualization,
- functional-group and SMILES-linked ligand annotations,
- PyMOL/script/session generation,
- PROTACability-style degrader-readiness triage, and
- downstream handoff to PROTAC Builder or related design workflows.

V-LiSEMOD is intended as a **low-maintenance lab handoff tool** for users who need to explore viral ligand structures without directly working inside the database or codebase.

---

<a id="why-vlisemod"></a>

## 🎯 Why V-LiSEMOD?

Viral protein–ligand co-crystal structures contain valuable information for antiviral drug discovery, inhibitor optimization, warhead design, and degrader-oriented medicinal chemistry.

However, the relevant evidence is often scattered across multiple structural and chemical layers:

```text
PDB structure
    ↓
ligand identity
    ↓
protein target / virus / chain context
    ↓
ligand atoms and functional groups
    ↓
protein-ligand contacts
    ↓
solvent exposure
    ↓
SMILES and atom mappings
    ↓
candidate linker atoms / warhead linkability
    ↓
target lysine accessibility
    ↓
PROTACability-style triage
```

V-LiSEMOD brings these layers into a searchable, exportable, and visually interpretable web platform.

The central design questions are:

```text
Which viral protein-ligand structures are worth reviewing?
Which ligand atoms are exposed enough to modify?
Which atoms appear important for binding and should be protected?
Which targets and structures may be better candidates for degrader-style design?
```

---

<a id="repository-navigation"></a>

## 🧭 Repository Navigation

<p align="center">
  <a href="#quick-start">
    <img src="https://img.shields.io/badge/Quick%20Start-Run%20Locally-orange?style=for-the-badge&logo=python" alt="Run locally">
  </a>
  <a href="#application-pages">
    <img src="https://img.shields.io/badge/Application%20Pages-User%20Guide-22c55e?style=for-the-badge" alt="Application pages">
  </a>
  <a href="#data-policy">
    <img src="https://img.shields.io/badge/Data%20Policy-Keep%20Repo%20Light-lightgrey?style=for-the-badge&logo=github" alt="Data policy">
  </a>
  <a href="#troubleshooting">
    <img src="https://img.shields.io/badge/Troubleshooting-Common%20Fixes-red?style=for-the-badge" alt="Troubleshooting">
  </a>
</p>

- [Overview](#overview)
- [Why V-LiSEMOD?](#why-vlisemod)
- [Core capabilities](#core-capabilities)
- [Application pages](#application-pages)
- [Conceptual workflow](#conceptual-workflow)
- [Scientific interpretation](#scientific-interpretation)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
- [Environment configuration](#environment-configuration)
- [Database and data layers](#database-and-data-layers)
- [PROTACability data regeneration](#protacability-data-regeneration)
- [Deployment](#deployment)
- [GitHub workflow](#github-workflow)
- [Data and output policy](#data-policy)
- [Troubleshooting](#troubleshooting)
- [Developer notes](#developer-notes)
- [Roadmap](#roadmap)
- [Citation](#citation)
- [Contact](#contact)

---

<a id="core-capabilities"></a>

## ✨ Core Capabilities

| Capability | Purpose |
|---|---|
| **Home structure explorer** | Select virus, PDB, ligand, and chain context for structure-specific review. |
| **Protein Query** | Search structures by virus, protein type, ligand, and dataset category. |
| **Ligand Indexer** | Start from a ligand code or synonym and locate mapped PDB-chain-residue contexts. |
| **Ligand Comparison** | Compare the same ligand across multiple structures using interaction and distance profiles. |
| **Ligand SVG / atom visualization** | Generate visual ligand outputs and solvent-exposed atom views. |
| **PyMOL/session export** | Create PyMOL-ready scripts and sessions for manual structural inspection. |
| **Interaction diagrams** | Review ligand-protein contact distributions and atom-level interaction burden. |
| **PROTACability Assessment** | Triage targets and structures using transparent degrader-readiness evidence layers. |
| **Warhead linkability evidence** | Evaluate solvent-exposed, chemically modifiable ligand atoms that may support linker attachment. |
| **Target lysine accessibility evidence** | Separately evaluate target-side lysine availability as a ubiquitination-readiness cue. |
| **Export workflows** | Download selected datasets as CSVs, Excel workbooks, or ZIP bundles depending on page. |
| **Optional Drug GPT / BioGPT module** | Disabled by default for lightweight deployment; can be re-enabled when model resources are ready. |

---

<a id="application-pages"></a>

## 🧩 Application Pages

### Home

The main structure-specific entry point.

**Typical use:** choose a virus, PDB code, ligand, and chain context, then generate visual or structural outputs.

**Capabilities:**

- select virus, PDB, ligand, and chain context,
- generate ligand SVGs,
- create solvent-exposed atom visualizations,
- generate PyMOL scripts/sessions,
- inspect structure-specific ligand context.

**Primary data dependencies:**

- `ligand_atoms`
- `Functional_Group_Atoms`
- `receptor_binding_pocket`
- `RUPLEY_SASA_DATA`
- `Ligand_Arp_Diagram`
- `Functional_GROUPED`

---

### About

The scientific motivation and project context page.

**Typical use:** explain V-LiSEMOD to collaborators, reviewers, or lab users.

**Capabilities:**

- describe the viral ligand and degrader-design motivation,
- connect the platform to the broader tool ecosystem,
- explain why solvent-exposed ligand atoms and warhead/linkability evidence matter.

**Primary data dependencies:**

- template-driven page,
- no heavy scientific table required.

---

### Protein Query

The main target-centric search and export workflow.

**Typical use:** filter viral structures by virus, protein type, ligand, or dataset category.

**Capabilities:**

- search by virus,
- filter by protein class/type,
- filter by ligand,
- export selected datasets,
- create CSV/Excel/ZIP bundles,
- open PROTACability summary popups,
- route into the full PROTACability dashboard.

**Primary data dependencies:**

- `Virus_Proteins`
- `ligand_synonyms`
- PROTACability tables when imported.

---

### Ligand Indexer

The ligand-first lookup page.

**Typical use:** start with a ligand code or synonym, then find PDB-chain-residue contexts where the ligand appears.

**Capabilities:**

- search by ligand code,
- search by ligand synonym,
- select a PDB-chain-residue instance,
- generate interaction chart snapshots,
- inspect ligand-specific structural contexts.

**Primary data dependencies:**

- `ligand_synonyms`
- `Ligand_Arp_Diagram`
- `Arpeggio_Contacts_Data`

---

### Ligand Comparison

The multi-structure ligand comparison workflow.

**Typical use:** compare one ligand across multiple PDB-chain-residue contexts.

**Capabilities:**

- compare interaction type distributions,
- inspect ligand distance profiles,
- review atom-level interaction burden,
- identify conserved versus context-specific interactions,
- reason about potentially modifiable ligand atoms.

**Primary data dependencies:**

- `Ligand_Atoms_Smiles`
- `SMILES_MAP_PDB`
- `Arpeggio_Contacts_Data`
- `ligand_synonyms`

---

### PROTACability Assessment

The degrader-readiness dashboard.

**Typical use:** evaluate target/structure/ligand combinations for transparent degrader-design triage.

**Capabilities:**

- browse Target Browser, Protein Summary, Structure Summary, and Chain Details,
- filter by score, tier, ligand context, evidence level, SMILES source, and atom evidence,
- inspect warhead linkability evidence,
- inspect target lysine accessibility evidence,
- review protein structural priority,
- review ternary geometry cues,
- export filtered or full PROTACability tables.

**Primary data dependencies:**

- `protacability_assessment`
- `protacability_lysine_proximity`
- `protacability_ligand_inventory`
- `protacability_warhead_linkability`
- `protacability_degrader_readiness`

---

### PROTAC Builder Handoff

A downstream modular degrader-design workflow.

**Typical use:** use V-LiSEMOD to identify promising viral warheads or ligand contexts, then move into a builder workflow.

**Capabilities:**

- link out to PROTAC Builder,
- use V-LiSEMOD warhead/linkability evidence as starting design rationale,
- explore linker/recruiter combinations downstream.

**Primary data dependencies:**

- external module or external URL,
- V-LiSEMOD provides conceptual and structure-guided handoff.

---

### Drug GPT / BioGPT

An optional assistant-style module.

**Typical use:** disabled by default for lightweight lab deployment.

**Capabilities when intentionally enabled:**

- route into a local assistant module,
- connect molecular/viral data to a language-model interface,
- support future trained-model workflows.

**Default behavior:**

This module should remain disabled unless model files, server memory, and deployment configuration are ready.

**Related endpoint note:**

If a template uses:

```html
{{ url_for('dp.home') }}
```

but the Drug GPT blueprint is disabled, Flask may raise:

```text
BuildError: Could not build url for endpoint 'dp.home'.
```

A safe disabled-state route is:

```html
{{ url_for('drug_gpt_disabled_home') }}
```

Use the real blueprint endpoint only when the `dp` blueprint is registered and active.

---

<a id="conceptual-workflow"></a>

## 🧬 Conceptual Workflow

A typical V-LiSEMOD workflow looks like this:

```text
1. Start from a virus, protein, PDB ID, or ligand
2. Filter relevant structures and ligand contexts
3. Inspect protein-ligand interactions
4. Review ligand atom exposure and functional-group evidence
5. Compare ligand behavior across structures
6. Evaluate PROTACability evidence layers
7. Export CSV / Excel / ZIP outputs
8. Hand off promising ligands or warheads to downstream design tools
```

A degrader-focused workflow looks like this:

```text
Viral protein-ligand co-crystal
        ↓
Ligand contact and solvent-exposure analysis
        ↓
Candidate modifiable ligand atoms
        ↓
Warhead linkability evidence
        ↓
Target lysine accessibility evidence
        ↓
Protein structural priority
        ↓
PROTACability triage
        ↓
PROTAC Builder / modeling / experimental design
```

---

<a id="scientific-interpretation"></a>

## 🧠 Scientific Interpretation

V-LiSEMOD provides **computational triage and structural decision support**, not experimental proof.

The PROTACability and degrader-readiness outputs are transparent, interpretable heuristics. They are intended to support hypothesis generation and prioritization.

They should **not** be interpreted as experimentally validated degradation predictions.

### Preferred interpretation language

Use wording like:

> **Warhead linkability evaluates whether the bound ligand contains solvent-exposed, chemically modifiable atoms that may tolerate linker attachment while preserving binding. Target lysine accessibility is evaluated separately as a ubiquitination-readiness cue. Ternary geometry cues are hypothesis-generating only and should not be interpreted as proof of productive degradation.**

### Avoid wording like

```text
lysine-linker site
linker attaches to target lysine
lysine linkability
this target is guaranteed to degrade
validated degrader prediction
```

### PROTACability evidence layers

| Evidence layer | Meaning |
|---|---|
| **Overall Degrader Readiness** | Combined triage score summarizing degrader-design evidence. |
| **Protein Structural Priority** | Protein/structure-level score derived from ligand context, structural context, lysine availability, pI, and related evidence. |
| **Warhead Linkability** | Ligand-centered evidence that a bound ligand may contain solvent-exposed, chemically modifiable atoms. |
| **Target Lysine Accessibility** | Target-side evidence that accessible lysines are present as ubiquitination-readiness cues. |
| **Ternary Geometry Cue** | Weak, hypothesis-only cue based on ligand-proximal exposed lysines. |

---

<a id="repository-layout"></a>

## 📦 Repository Layout

The exact repository structure may evolve, but a typical V-LiSEMOD layout is:

```text
VLISEMOD/
├── app.py                              # Main Flask application entry point
├── DRUGapp.py                          # Optional Drug GPT / BioGPT module
├── requirements.txt                    # Python dependencies
├── templates/                          # Flask templates
│   ├── base.html
│   ├── index.html
│   ├── about.html
│   ├── protein_query.html
│   ├── ligand_indexer.html
│   ├── compare_ligands.html
│   └── protacability_page.html
├── static/                             # CSS, JavaScript, images, generated frontend assets
│   ├── css/
│   ├── js/
│   ├── charts/
│   ├── coordinate_cache/
│   ├── ligand_images/
│   └── ligand_sdf_cache/
├── TOOLS/                              # Importers, maintenance utilities, and data scripts
│   └── import_protacability_data.py
├── PDB_FILES/                          # Large structure files; keep out of Git
├── output_csvs/                        # Generated exports; keep out of Git
├── pml_sessions/                       # Generated PyMOL sessions; keep out of Git
├── viral_data.db                       # Primary SQLite database; usually not committed
├── users.db                            # Optional user/auth database; do not commit if private
└── README.md                           # Project documentation
```

---

<a id="quick-start"></a>

## ⚡ Quick Start

Clone the repository:

```bash
git clone <YOUR_REPOSITORY_URL>
cd VLISEMOD
```

Create and activate an environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

Or use conda:

```bash
conda create -n viraldb python=3.10
conda activate viraldb
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Create a local environment file if needed:

```bash
cp .env.example .env
```

Run locally:

```bash
python app.py
```

Or with Waitress:

```bash
waitress-serve --listen=127.0.0.1:5002 app:app
```

Open:

```text
http://127.0.0.1:5002
```

If `python app.py` uses a different port, use the URL shown in the terminal.

---

<a id="environment-configuration"></a>

## 🔐 Environment Configuration

Recommended environment variables:

| Variable | Purpose |
|---|---|
| `FLASK_SECRET_KEY` | Flask session secret. Required for production-style deployment. |
| `ENABLE_DRUG_GPT` | Enables Drug GPT routes when set intentionally. |
| `ENABLE_LOCAL_LLM` | Allows local model loading. Keep off unless model resources are ready. |
| `SHOW_DRUG_GPT_NAV` | Shows or hides Drug GPT from navigation. |
| `PROTAC_BUILDER_EXTERNAL_URL` | Overrides the external PROTAC Builder URL. |
| `MAX_SESSIONS_PER_USER` | Controls per-user session concurrency if enabled. |
| `SESSION_CONCURRENCY_MODE` | Controls behavior when session limits are reached. |
| `DATABASE_URL` | Optional external database connection string for deployed environments. |
| `PORT` | Runtime port used by deployment platforms. |

Example `.env` starter:

```bash
FLASK_SECRET_KEY=replace-me-with-a-real-secret
ENABLE_DRUG_GPT=false
ENABLE_LOCAL_LLM=false
SHOW_DRUG_GPT_NAV=false
PROTAC_BUILDER_EXTERNAL_URL=https://protacbuilder.com
MAX_SESSIONS_PER_USER=3
SESSION_CONCURRENCY_MODE=warn
PORT=5002
```

For lightweight deployment, keep:

```bash
ENABLE_DRUG_GPT=false
ENABLE_LOCAL_LLM=false
SHOW_DRUG_GPT_NAV=false
```

This prevents the app from attempting to load or download large language-model checkpoints.

---

<a id="database-and-data-layers"></a>

## 🧾 Database and Data Layers

The primary application database is typically:

```text
viral_data.db
```

Important tables/data layers include:

| Table / data layer | Purpose |
|---|---|
| `Virus_Proteins` | Virus and protein target metadata. |
| `ligand_synonyms` | Ligand code and synonym lookup. |
| `ligand_atoms` | Ligand atom records and structure-linked atom data. |
| `Ligand_Atoms_Smiles` | Ligand SMILES-linked atom information. |
| `SMILES_MAP_PDB` | PDB-to-SMILES atom mapping. |
| `Arpeggio_Contacts_Data` | Protein-ligand contact data. |
| `Ligand_Arp_Diagram` | Ligand interaction diagram data. |
| `RUPLEY_SASA_DATA` | Shrake-Rupley-style solvent accessibility records. |
| `Functional_Group_Atoms` | Functional-group atom annotations. |
| `Functional_GROUPED` | Grouped functional annotation data. |
| `receptor_binding_pocket` | Binding-pocket residue context. |
| `protacability_assessment` | Core PROTACability assessment table. |
| `protacability_lysine_proximity` | Ligand-proximal lysine and accessibility cues. |
| `protacability_ligand_inventory` | Ligand inventory for PROTACability workflows. |
| `protacability_warhead_linkability` | Warhead/linkability enrichment evidence. |
| `protacability_degrader_readiness` | Combined degrader-readiness output table. |

---

<a id="protacability-data-regeneration"></a>

## 🧪 PROTACability Data Regeneration

When PROTACability tables need to be regenerated, run the workflow scripts in order.

```bash
python 00_Protein_Expansion.py
python 00_Protein_Expansion2.py
python 00_Protein_Expansion3.py --fresh --workers 8
python 01_PROTACability_Warhead_Linkability_Enrichment.py --workers 8 --component-smiles Components-smiles-stereo-oe.smi
python TOOLS/import_protacability_data.py
```

Then verify the app:

```text
/
 /about
 /query_protein_virus_page
 /ligand_indexer
 /compare_ligands
 /protacability_page
 /healthz
```

Checklist:

- homepage renders,
- navigation links work,
- Protein Query exports download,
- Ligand Indexer searches work,
- Ligand Comparison charts render,
- PROTACability filters populate,
- PROTACability detail modal opens,
- Drug GPT remains disabled unless intentionally enabled,
- `/healthz` returns successfully.

---

<a id="deployment"></a>

## 🚀 Deployment

V-LiSEMOD can be deployed from GitHub to Heroku or another Flask-compatible platform.

### Recommended GitHub-first deployment flow

Instead of pushing directly to Heroku Git, push to GitHub first:

```bash
git status
git add .
git commit -m "Update V-LiSEMOD app"
git push origin main
```

Then connect the GitHub repository inside the Heroku dashboard:

```text
Heroku app
  → Deploy
  → Deployment method: GitHub
  → Connect repository
  → Select branch
  → Deploy Branch
```

Automatic deploys can be enabled after the manual deployment works.

### Heroku CLI note

If this appears locally:

```bash
heroku: command not found
```

that does not mean the Heroku account is wrong. It means the Heroku CLI is not installed or not available in the active terminal PATH.

You can still deploy through the Heroku dashboard using GitHub integration.

### Production process command

A typical `Procfile` for Flask + Waitress may look like:

```Procfile
web: waitress-serve --port=$PORT app:app
```

If using Gunicorn instead:

```Procfile
web: gunicorn app:app
```

Use whichever server is present in `requirements.txt` and validated for the deployment target.

---

<a id="github-workflow"></a>

## 🔁 GitHub Workflow

Recommended development loop:

```bash
git status
python app.py
```

After testing locally:

```bash
git add .
git commit -m "Describe the change"
git push origin main
```

If the remote branch has newer commits:

```bash
git pull --rebase origin main
git push origin main
```

Before pushing, check what will be committed:

```bash
git status --short
git diff --stat
```

Avoid committing large generated files:

```bash
find . -type f -size +50M \
  -not -path "./.git/*" \
  -exec ls -lh {} \;
```

---

<a id="data-policy"></a>

## 🧹 Data and Output Policy

This repository should remain lightweight and source-code focused.

### Recommended to commit

- Flask application source code,
- templates,
- static frontend assets,
- lightweight JavaScript and CSS,
- route helpers,
- import scripts,
- documentation,
- `.env.example`,
- small example files if explicitly intended.

### Do not commit

```text
PDB_FILES/
viral_data.db
users.db
Components-smiles-stereo-oe.smi
output_csvs/
pml_sessions/
static/coordinate_cache/
static/ligand_sdf_cache/
static/charts/
static/ligand_images/ when generated locally
.env
model files
model checkpoints
API keys
Heroku tokens
private credentials
large generated ZIP files
```

### Useful checks

Check repository size:

```bash
du -sh .
du -h --max-depth=1 . | sort -hr
```

Find large files:

```bash
find . -type f -size +50M \
  -not -path "./.git/*" \
  -exec ls -lh {} \;
```

Confirm ignored files are ignored:

```bash
git check-ignore -v viral_data.db
git check-ignore -v .env
git check-ignore -v PDB_FILES/example.pdb
git check-ignore -v output_csvs/example.csv
```

---

<a id="troubleshooting"></a>

## 🧯 Troubleshooting

### Homepage crashes with `BuildError: Could not build url for endpoint 'dp.home'`

Cause:

The template is trying to build a link to the Drug GPT blueprint:

```html
{{ url_for('dp.home') }}
```

but that blueprint is not currently registered.

Safe disabled-state fix:

```html
{{ url_for('drug_gpt_disabled_home') }}
```

Use `dp.home` only when the Drug GPT blueprint is active and registered.

To inspect active Flask routes:

```bash
python -m flask --app app routes
```

or:

```bash
flask --app app routes
```

---

### App starts slowly or tries to download a huge model

Cause:

Drug GPT or local LLM loading is probably enabled.

Fix:

```bash
ENABLE_DRUG_GPT=false
ENABLE_LOCAL_LLM=false
SHOW_DRUG_GPT_NAV=false
```

Then restart the app.

---

### `heroku: command not found`

Cause:

Heroku CLI is not installed or not on `PATH`.

Options:

1. Deploy through GitHub integration in the Heroku dashboard.
2. Install the Heroku CLI if command-line Heroku operations are needed.

The GitHub-first workflow is usually cleaner for this project.

---

### PROTACability filters are empty

Cause:

PROTACability tables may not have been imported.

Fix:

```bash
python TOOLS/import_protacability_data.py
```

Then restart the app.

---

### Ligand images or charts are missing

Possible causes:

- generated static caches are absent,
- chart output folders were cleared,
- structure or ligand-specific output was never generated.

Fix:

- regenerate the relevant page output,
- clear stale caches if needed,
- confirm write permissions for `static/charts/`, `static/ligand_images/`, and related generated folders.

---

### Navigation links fail after route changes

Run:

```bash
python -m flask --app app routes
```

Then compare active endpoint names to all `url_for(...)` calls in templates.

Search templates:

```bash
grep -R "url_for" templates
```

---

<a id="developer-notes"></a>

## 🛠️ Developer Notes

### Design goals

V-LiSEMOD should be:

- easy for lab members to use,
- safe to deploy without local LLM/model dependencies,
- organized around clear biological workflows,
- able to export useful datasets,
- modular enough to connect with downstream tools,
- interpretable rather than black-box.

### Optional modules

Optional modules like Drug GPT / BioGPT should remain disabled unless:

- the model files are available,
- deployment memory is sufficient,
- routes are registered intentionally,
- environment variables are set correctly,
- the navigation and templates are updated accordingly.

### Route sanity check

After any template or navigation change:

```bash
python -m flask --app app routes
```

Then visit:

```text
/
 /about
 /query_protein_virus_page
 /ligand_indexer
 /compare_ligands
 /protacability_page
 /healthz
```

### Recommended pre-commit test

```bash
python app.py
```

Open the homepage and click through the navigation before pushing.

---

<a id="roadmap"></a>

## 🧭 Roadmap

Potential development directions:

- stronger route/API validation tests,
- GitHub Actions smoke tests,
- Docker-based deployment,
- clearer example dataset bundle,
- manuscript-quality report exports,
- improved candidate linker atom visualization,
- 3D candidate atom highlighting in the viewer,
- better PROTACability detail report pages,
- downloadable per-structure evidence reports,
- tighter PROTAC Builder integration,
- improved dark/light mode consistency across all pages,
- polished graph scaling for multi-selection ligand comparison,
- optional authentication for private lab deployments,
- reintroduction of a trained, deployment-safe assistant module after validation.

---

<a id="manuscript-positioning"></a>

## 📄 Manuscript Positioning

Working title:

> **V-LiSEMOD: A Viral Ligand Solvent-Exposed Moiety Database for Structure-Guided Warhead Discovery and Degrader-Readiness Triage**

Potential manuscript framing:

```text
V-LiSEMOD organizes viral protein-ligand co-crystal structures into searchable,
visual, and exportable workflows that support structure-guided inhibitor analysis,
warhead discovery, solvent-exposed ligand atom identification, and transparent
PROTACability-style degrader-readiness triage.
```

Suggested feature-demonstration sections:

1. V-LiSEMOD organizes viral protein-ligand structures into searchable workflows.
2. Ligand-centered visualization highlights solvent-exposed atoms and local structural context.
3. Ligand Comparison identifies conserved and context-specific interaction burdens.
4. PROTACability separates warhead evidence from target lysine accessibility.
5. PROTACability identifies representative structures for downstream design.
6. V-LiSEMOD supports downstream PROTAC Builder handoff.

---

<a id="citation"></a>

## 🧬 Citation

A formal manuscript or software citation can be added here when available.

For now, cite the repository and platform:

```text
Schulz, J.-M. V-LiSEMOD: Viral Ligand Solvent-Exposed Moiety Database for Structure-Guided Warhead Discovery and Degrader-Readiness Triage. GitHub repository: <repository URL>. Web platform: <deployment URL>.
```

If citing the scientific motivation for degrader-readiness and targeted degrader thinking, include the related targeted degrader publication where appropriate:

```text
Khurshid et al. Targeted degrader technologies as prospective SARS-CoV-2 therapies. Drug Discovery Today. 2024.
```

---

<a id="contact"></a>

## 📬 Contact

For questions, workflow support, bug reports, or collaboration inquiries:

<p align="center">
  <a href="mailto:jxs794@miami.edu?subject=V-LiSEMOD%20Question%20%2F%20Collaboration">
    <img src="https://img.shields.io/badge/Joseph--Michael%20Schulz-jxs794%40miami.edu-blue?style=for-the-badge&logo=gmail" alt="Email Joseph-Michael Schulz">
  </a>
</p>

---

<a id="repository-description"></a>

## 🧾 Repository Description

> V-LiSEMOD is a Flask-based viral ligand structure exploration platform for protein-ligand interaction analysis, solvent-exposed ligand atom discovery, warhead linkability evidence review, and PROTACability-style degrader-readiness triage.

---

<a id="practical-takeaway"></a>

## 🙌 Practical Takeaway

Use V-LiSEMOD when you need to move from:

```text
viral protein-ligand structure
```

to:

```text
interaction evidence + solvent-exposed ligand atoms + candidate warhead/linker rationale + PROTACability triage + downstream design handoff
```

The platform helps convert viral structural data into interpretable medicinal chemistry and degrader-design evidence while remaining transparent about its heuristic, hypothesis-generating nature.
