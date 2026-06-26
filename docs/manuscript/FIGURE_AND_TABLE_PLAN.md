# Figure And Table Plan

## Figures

## Figure-package status update

- Public screenshot evidence now exists for Warhead Hunter, V-LiSEMOD, PROTAC Builder, E3 Ligandalyzer, PyMACS, and Warhead Hunter API docs.
- A rough ecosystem composite exists at `docs/manuscript/figures/annotated/fig01_ecosystem_overview_draft.png`.
- Cropped manuscript-oriented panels now exist for Warhead Hunter, V-LiSEMOD, PROTAC Builder, E3 Ligandalyzer, and PyMACS.
- The AutoDock-Vina PrepServer figure slot remains pending because the current workspace is V-LiSEMOD rather than the local AutoDock application named in the prompt.
- Several long raw screenshots need crop-based use because sticky or dynamic layouts repeated sections during full-page capture.
- E3 Ligandalyzer and PROTAC Builder captures are usable for draft figures but may merit recapture for cleaner no-overlay or post-modal panels.

### Figure 1. Platform overview schematic

- Purpose: show the end-to-end logic from viral protein-ligand structure input to evidence review and companion-tool handoff.
- Source material needed: application modules, data layers, and conceptual workflow from current docs.
- Evidence already exists: yes, at the workflow-description level.
- Screenshot or data needed: one clean schematic rather than raw screenshots.
- Draft caption: “Overview of the V-LiSEMOD workflow from curated viral co-crystal context through ligand-centered evidence review, solvent-exposed atom interpretation, PROTACability-style triage, and downstream companion-tool handoff.”
- Overclaiming risk: moderate if the schematic is drawn like an automated prediction pipeline rather than a guided review workflow.

### Figure 2. Application map

- Purpose: orient readers to the main pages and user journeys.
- Source material needed: home page, Protein Query, Ligand Indexer, Ligand Comparison, PROTACability Assessment, optional Drug GPT route.
- Evidence already exists: yes.
- Screenshot or data needed: representative screenshots or a composite panel.
- Draft caption: “Map of the main V-LiSEMOD modules for structure-specific, target-centric, ligand-centric, and heuristic triage workflows.”
- Overclaiming risk: low if optional modules are clearly marked optional.

### Figure 3. Data architecture figure

- Purpose: explain how Flask routes, `viral_data.db`, generated caches, templates, and optional RANDY mode fit together.
- Source material needed: `app.py`, `RANDY/`, docs, and generated-folder inventory.
- Evidence already exists: yes.
- Screenshot or data needed: architecture schematic.
- Draft caption: “Simplified V-LiSEMOD architecture showing local SQLite-backed operation, generated runtime assets, and optional RANDY-backed route groups.”
- Overclaiming risk: medium if drawn like a hardened production architecture.

### Figure 4. Ligand-centered evidence figure

- Purpose: show how ligand atoms, interactions, SASA, and functional-group annotations converge in one workflow.
- Source material needed: structure explorer outputs and ligand imagery/chart outputs.
- Evidence already exists: partially.
- Screenshot or data needed: one curated example structure with exposed atoms and interaction overlays.
- Draft caption: “Example ligand-centered evidence review combining atom mapping, interaction context, solvent exposure, and functional-group annotations.”
- Overclaiming risk: low if framed as an example workflow rather than universal performance.

### Figure 5. PROTACability interpretation figure

- Purpose: explain warhead linkability, lysine accessibility, structural priority, and combined readiness heuristics.
- Source material needed: PROTACability page plus imported CSV/database evidence layers.
- Evidence already exists: yes, at the software and table level.
- Screenshot or data needed: a representative dashboard view or a schematic evidence stack.
- Draft caption: “Interpretation framework for V-LiSEMOD PROTACability outputs, emphasizing structural-priority triage rather than experimental degradation validation.”
- Overclaiming risk: high unless the heuristic warning is explicit in the caption and figure text.

### Figure 6. Companion-tool ecosystem figure

- Purpose: place V-LiSEMOD within the broader induced-proximity and structure-guided design ecosystem.
- Source material needed: README positioning and external companion-tool references.
- Evidence already exists: yes, conceptually.
- Screenshot or data needed: simple ecosystem diagram.
- Draft caption: “Positioning of V-LiSEMOD alongside companion tools supporting downstream warhead, PROTAC, and E3-ligand exploration workflows.”
- Overclaiming risk: medium if the figure implies deeper technical integration than currently documented.

### Figure 7. Public API and reproducibility evidence figure

- Purpose: document that selected public documentation and lightweight endpoint surfaces are reachable.
- Source material needed: Warhead Hunter API docs, small public endpoint captures, and V-LiSEMOD `/healthz`.
- Evidence already exists: yes.
- Screenshot or data needed: one API-doc screenshot plus small endpoint evidence files.
- Draft caption: “Public documentation and lightweight endpoint checks supporting software-availability and reproducibility-oriented manuscript evidence.”
- Overclaiming risk: medium if framed as production-scale API validation rather than lightweight public evidence.

## Tables

### Table 1. Platform capabilities table

- Purpose: summarize the main modules and what each module supports.
- Source material needed: `docs/APP_GUIDE.md`.
- Evidence already exists: yes.
- Screenshot or data needed: none.
- Draft caption: “Current V-LiSEMOD modules and their primary user-facing roles.”
- Overclaiming risk: low.

### Table 2. Data-layer summary

- Purpose: list the main database tables and their analytical roles.
- Source material needed: `docs/DATABASE.md` and inspected SQLite schema.
- Evidence already exists: yes.
- Screenshot or data needed: none.
- Draft caption: “Major V-LiSEMOD data layers supporting structure, ligand, interaction, solvent-exposure, and PROTACability workflows.”
- Overclaiming risk: low.

### Table 3. PROTACability evidence-layer interpretation table

- Purpose: distinguish each heuristic evidence layer from experimental claims.
- Source material needed: `docs/PROTACABILITY.md` and manuscript claims matrix.
- Evidence already exists: yes.
- Screenshot or data needed: none.
- Draft caption: “Interpretation guardrails for V-LiSEMOD PROTACability-style evidence layers.”
- Overclaiming risk: high if the table is not explicit about what the outputs do not prove.

### Table 4. App module and user-question table

- Purpose: connect modules to the practical questions they help answer.
- Source material needed: `docs/APP_GUIDE.md`.
- Evidence already exists: yes.
- Screenshot or data needed: none.
- Draft caption: “Representative research questions supported by each V-LiSEMOD module.”
- Overclaiming risk: low.

### Table 5. Limitations and future-work table

- Purpose: keep limitations visible and manuscript-safe.
- Source material needed: claims matrix and reproducibility plan.
- Evidence already exists: yes.
- Screenshot or data needed: none.
- Draft caption: “Current limitations and future-work directions for V-LiSEMOD.”
- Overclaiming risk: low if limitations remain specific.

### Table 6. Reproducibility checklist table

- Purpose: summarize what a reviewer or collaborator needs to reproduce a run.
- Source material needed: deployment guide and validation plan.
- Evidence already exists: yes.
- Screenshot or data needed: none.
- Draft caption: “Checklist for reproducible local or provisioned V-LiSEMOD runs.”
- Overclaiming risk: medium if the table implies turnkey reproducibility without data provisioning.
