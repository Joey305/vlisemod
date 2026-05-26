# V-LiSEMOD Site Documentation

## 1. Executive Summary
V-LiSEMOD is a viral ligand solvent-exposed moiety database and structure-guided analysis platform for exploring curated viral protein-ligand structures, ligand interactions, solvent-exposed atoms, ligand-centered warhead/linkability reasoning, and PROTAC/degrader-readiness triage. The application combines a Flask web app, a SQLite database (`viral_data.db`), structure-derived caches and generated exports, and a set of enrichment/import scripts that support both direct user workflows and downstream design handoff.

At a practical level, V-LiSEMOD helps a user:

- Start from a virus, protein class, or ligand.
- Inspect experimentally resolved protein-ligand contexts.
- Export PyMOL sessions, CSV tables, and ligand-centered visualization assets.
- Compare interaction patterns across structures.
- Review PROTACability-style structural-priority evidence for viral targets.
- Hand promising warhead candidates downstream into the standalone PROTAC Builder workflow.

Important interpretation note: PROTACability outputs in V-LiSEMOD are transparent structural-priority and design-readiness heuristics. They are useful for hypothesis generation and triage, but they are not experimentally validated degradation predictions and should not be interpreted as proof of productive targeted degradation.

## 2. Table of Contents
- [1. Executive Summary](#1-executive-summary)
- [2. Table of Contents](#2-table-of-contents)
- [3. High-Level Application Map](#3-high-level-application-map)
- [4. User Workflow Overview](#4-user-workflow-overview)
- [5. Page-by-Page Documentation](#5-page-by-page-documentation)
  - [Home](#home)
  - [About](#about)
  - [Protein Query](#protein-query)
  - [Ligand Indexer](#ligand-indexer)
  - [Ligand Comparison](#ligand-comparison)
  - [PROTACability Assessment](#protacability-assessment)
  - [PROTAC Builder](#protac-builder)
  - [Drug GPT / BioGPT](#drug-gpt--biogpt)
- [6. Database Overview](#6-database-overview)
- [7. File/Folder Structure](#7-filefolder-structure)
- [8. Deployment / Environment Notes](#8-deployment--environment-notes)
- [9. Maintenance Checklist](#9-maintenance-checklist)
- [10. Known Limitations](#10-known-limitations)
- [11. Future Work](#11-future-work)
- [12. Quick Start for a New Lab User](#12-quick-start-for-a-new-lab-user)
- [13. Quick Start for a Maintainer](#13-quick-start-for-a-maintainer)

## 3. High-Level Application Map
| Page / Module | URL / Route | Primary Purpose | Main Inputs | Main Outputs | Database Tables / Files Used | Maintenance Notes |
|---|---|---|---|---|---|---|
| Home | `/` | Landing page, main session builder, workflow overview | Virus, PDB code, ligand, optional feature checkboxes | Downloadable PyMOL session, ligand image generation handoff | `ligand_atoms`, `Functional_Group_Atoms`, `receptor_binding_pocket`, `distal_atoms`, `RUPLEY_SASA_DATA`, `Ligand_Arp_Diagram`, `Functional_GROUPED`; generated `pml_sessions/`, `static/ligand_images/` | Main “session builder” workflow. Template still contains a Drug GPT card, so keep hidden-module copy in sync with deployment flags. |
| About | `/about` | Scientific framing and project context | None | Informational content and external links | Template-driven only | Keep scientific framing current with manuscript/review language. |
| Protein Query | `/query_protein_virus_page` | Filter by virus, protein, ligand, and export datasets | Virus names, protein types, optional ligand, dataset checkboxes | ZIP containing CSVs plus combined Excel workbook | `Virus_Proteins`, `ligand_synonyms`, export dataset tables from `data_set_queries`; PROTACability tables if present | Main route for export-ready filtered structure sets. |
| Ligand Display | returned by `/generate_ligand_images` | 2D/3D ligand display, SASA image toggle, PoseView diagram launch | Virus, PDB, ligand, chain/residue context | SVG downloads, 3D viewer, PoseView diagram, warhead handoff | `Functional_GROUPED`, `Ligand_Arp_Diagram`, `RUPLEY_SASA_DATA`, `SMILES_MAP_PDB`; generated `static/ligand_images/`, `static/coordinate_cache/`, `static/ligand_sdf_cache/` | Secondary workflow page, not in nav, but central to ligand visualization. |
| Ligand Indexer | `/ligand_indexer` | Ligand-first lookup and interaction chart generation | Ligand or synonym, PDB-chain-residue selection | Carousel of interaction charts | `ligand_synonyms`, `Ligand_Arp_Diagram`, `Arpeggio_Contacts_Data` | Older workflow retained and re-skinned into current shell. |
| Ligand Comparison | `/compare_ligands` | Multi-PDB ligand interaction comparison | Ligand or synonym, multiple PDB-chain-residue selections | Plotly charts, SVG view, PROTAC Builder handoff | `ligand_synonyms`, `Ligand_Atoms_Smiles`, `SMILES_MAP_PDB`, `Arpeggio_Contacts_Data` | Good downstream comparative page for warhead behavior and interaction burden. |
| PROTACability Assessment | `/protacability_page` | Multi-view degrader-readiness triage dashboard | Filters on virus, protein, readiness/warhead tiers, scores, ligand context, evidence flags | Paginated summaries, detail modal, 3D viewer, CSV export | `protacability_assessment`, `protacability_lysine_proximity`, `protacability_ligand_inventory`, optional `protacability_warhead_linkability`, optional `protacability_degrader_readiness`; `PDB_FILES/*.csv`; cached coordinates and ligand SDF | Most important maintenance surface for structural triage and enrichment integration. |
| PROTAC Builder handoff | External links via nav/footer/buttons; legacy redirect `/copy` | Downstream degrader design workflow | Selected ligand/warhead context from V-LiSEMOD | External builder launch | External site `https://protacbuilder.com` or `PROTAC_BUILDER_EXTERNAL_URL` | External and conceptually connected, not maintained inside this repo. |
| Tutorial banner | Included on Home, Protein Query, Ligand Comparison, Ligand Display | Lightweight help/tutorial link | Page-specific YouTube URL | Dismissible flyover banner | Template-only partial `_tutorial_flyover.html` | Page-local dismissal only; no persistence. |
| Footer / toolkit ecosystem | Shared in `base.html` | Cross-tool ecosystem framing | None | External links to Warhead Hunter, E3 Ligandalyzer, PROTAC Builder | Template-only | Keep links and descriptions synchronized with broader toolkit. |
| Drug GPT / BioGPT | `/drugapp/` when enabled, otherwise disabled stub | Hidden assistant / local LLM module | Free-text molecule question | Chat-style model response plus optional PubChem enrichment | `DRUGapp.py`, local model config, Hugging Face model loading when enabled | Hidden/disabled by default. Re-enable only intentionally. |
| Login/auth scaffold | `/login`, `/logout`, `/force-logout`, `/healthz` | User/session management | Email/password and session cookies | Auth/session records | `users.db` tables `users`, `sessions` | Auth exists, but current `EXEMPT_PATH_PREFIXES` includes `/`, which effectively exempts all routes from login checks. Verify before relying on auth as a gate. |

## 4. User Workflow Overview
Recommended day-to-day workflow:

1. Start at Home or About to understand the platform framing and choose a workflow.
2. Use Protein Query to filter viral structures by virus, protein class, ligand, and available datasets.
3. Use Ligand Indexer to inspect a ligand across structures and generate interaction chart snapshots.
4. Use Ligand Comparison to compare contacts, distances, and atom-level interaction burdens across selected PDB structures.
5. Use PROTACability Assessment to triage targets using protein structural priority, warhead linkability, target lysine accessibility, ternary geometry cue, and overall degrader-design readiness.
6. Use PROTAC Builder as the downstream modular design tool when a promising warhead/ligand context is identified.
7. Export CSVs, combined Excel bundles, ligand SVGs, SASA-highlighted SVGs, and PyMOL session outputs where available.

Interpretation guidance:

- Treat structure-derived outputs as experimentally anchored but still interpretive.
- Use solvent-exposed atoms, interaction preservation, and mapping coverage as support for chemistry reasoning, not as a guarantee of tractable optimization.
- Treat PROTACability/readiness outputs as transparent evidence layers rather than black-box predictions.

## 5. Page-by-Page Documentation

### Home
- Route: `/`
- Template: `templates/index.html`
- Purpose: Landing page, quick-start orientation, and main PyMOL session builder.
- What the user can do:
  - See a “Quick Start” modal with workflow choices.
  - Launch the main session builder for virus → PDB → ligand selection.
  - Generate a PyMOL session.
  - Trigger ligand image generation.
  - Jump to Protein Query, Ligand Comparison, upcoming features, or PROTAC Builder.
- Main controls/forms:
  - Virus dropdown populated from `/get_viruses`.
  - PDB dropdown populated from `/get_pdb_codes/<virus_name>`.
  - Ligand dropdown populated from `/get_ligands/<pdb_code>`.
  - Chain selector shown when multiple chains are present.
  - Optional checkboxes:
    - Binding Pocket
    - Shrake-Rupley SASA
    - Functional Group Atoms
- Backend routes/APIs called:
  - `/get_viruses`
  - `/get_pdb_codes/<virus_name>`
  - `/get_ligands/<pdb_code>`
  - `/check_functional_groups/<pdb_code>`
  - `/generate_pymol_session`
  - `/generate_ligand_images`
- Database tables used:
  - `ligand_atoms`
  - `Functional_Group_Atoms`
  - `receptor_binding_pocket`
  - `distal_atoms`
  - `RUPLEY_SASA_DATA`
  - `Ligand_Arp_Diagram`
  - `Functional_GROUPED`
- Output/export behavior:
  - Downloads a generated PyMOL script/session payload from `/generate_pymol_session`.
  - Opens a ligand image/3D inspection page when image generation succeeds.
- Scientific interpretation:
  - The home page is the entry point for structure-specific review, especially solvent exposure and local binding context.
- Known limitations:
  - Some chain-specific ligand mappings may be incomplete.
  - Functional group checkbox availability depends on `Functional_Group_Atoms` coverage.
- Maintenance notes:
  - Client behavior is driven by `static/js/scripts.js`.
  - The home template still contains a Drug GPT card linking to `dp.home`; if Drug GPT is fully disabled in a deployment, this card should stay hidden or be feature-flagged just like the nav.

### About
- Route: `/about`
- Template: `templates/about.html`
- Purpose: Explain the scientific motivation, project origin, and connected tool ecosystem.
- What the user can do:
  - Read the platform rationale.
  - See the link to the Drug Discovery Today review.
  - Discover connected tools such as Warhead Hunter, E3 Ligandalyzer, and PROTAC Builder.
- Main controls/forms:
  - None beyond navigation and external links.
- Backend routes/APIs called:
  - None beyond standard page render.
- Database tables used:
  - None directly.
- Output/export behavior:
  - None.
- Scientific interpretation:
  - The page explicitly frames V-LiSEMOD as a structure-first workspace for moving from resolved viral protein-ligand complexes toward warhead and degrader-oriented reasoning.
- Known limitations:
  - Narrative copy should be updated if project framing changes substantially.
- Maintenance notes:
  - Includes scientific origin language referencing the Drug Discovery Today paper. Keep citation wording aligned with lab preference and publication status.

### Protein Query
- Route: `/query_protein_virus_page`
- Template: `templates/query_protein_virus.html`
- Purpose: Filter viral structures by virus, protein class, optional ligand, and export-ready dataset selection.
- What the user can do:
  - Select one or more viruses.
  - Select one or more protein types.
  - Optionally constrain by ligand or ligand synonym.
  - See matching PDB choices.
  - Export multiple data products in one ZIP bundle.
  - Open the PROTACability popup and navigate into the full PROTACability dashboard.
- Main controls/forms:
  - Multi-select virus name.
  - Multi-select protein type.
  - Optional ligand filter.
  - Dataset export checkboxes:
    - Solvent Exposed Atoms
    - Ligand Atoms
    - Binding Pocket
    - Smiles and Functional Groups
    - Interatomic Interactions
    - Functional Group Atoms
    - Smiles & PDB Mapping
    - PROTACability Assessment
    - PROTACability Lysine Proximity
    - PROTACability Ligand Inventory
    - PROTACability Warhead Linkability when imported
    - PROTACability Degrader Readiness when imported
  - “Assess PROTACability” popup/button.
- Backend routes/APIs called:
  - `/get_virus_names_list_distinct`
  - `/get_protein_types_list_distinct`
  - `/get_ligands_with_synonyms`
  - `/get_pdbs_for_virus_protein`
  - `/export_data_to_excel`
  - `/protacability_page`
- Database tables used:
  - `Virus_Proteins`
  - `ligand_synonyms`
  - Export dataset source tables listed in `data_set_queries`
  - Optional PROTACability tables depending on import status
- Output/export behavior:
  - `/export_data_to_excel` builds CSV files for each selected dataset, writes a combined Excel workbook, packages both into `data_sets.zip`, then deletes the temporary output folder.
- Scientific interpretation:
  - This page is the best starting point for target-centric dataset export and structure bundle assembly.
  - The PROTACability popup is educational/navigation-focused and routes users into the dedicated dashboard rather than trying to explain all evidence layers inline.
- Known limitations:
  - Export coverage depends on source table completeness.
  - Ligand synonyms rely on the `ligand_synonyms` table.
  - Large structure selections may produce large ZIP bundles.
- Maintenance notes:
  - The PROTACability popup has explicit theme-aware styling for both light and dark mode and should stay aligned with shared CSS variables.
  - Tutorial/help banner comes from `_tutorial_flyover.html`.

### Ligand Indexer
- Route: `/ligand_indexer`
- Template: `templates/ligand_query.html`
- Purpose: Ligand-first lookup for mapped residues and interaction chart generation.
- What the user can do:
  - Search for a ligand by 3-letter code or synonym.
  - Select a PDB/chain/residue context where that ligand occurs.
  - Generate interaction visualizations in a carousel.
  - Open an interaction-definitions modal.
- Main controls/forms:
  - Ligand/synonym dropdown.
  - PDB-chain-residue dropdown.
  - “Show Interactions” button.
  - Definitions modal.
  - Chart carousel with save button.
- Backend routes/APIs called:
  - `/get_ligands_with_synonyms`
  - `/get_pdb_residue_by_ligand/<ligand_code>`
  - `/generate_charts`
  - `/drugapp/` HEAD check for chatbot availability
- Database tables used:
  - `ligand_synonyms`
  - `Ligand_Arp_Diagram`
  - `Arpeggio_Contacts_Data`
- Output/export behavior:
  - Returns chart image paths and displays them in a modal carousel.
  - Allows the currently visible chart image to be saved locally.
- Scientific interpretation:
  - Useful for quickly surveying how a ligand behaves across one structure-specific context at a time.
- Known limitations:
  - This page emphasizes chart output more than modern 3D context.
  - The hidden/disabled Drug GPT chatbot stub is still embedded in the UI shell.
- Maintenance notes:
  - This is an older workflow retained in the current design language; if simplifying the app in the future, compare its role with Ligand Comparison before removing it.

### Ligand Comparison
- Route: `/compare_ligands`
- Template: `templates/compare_ligands.html`
- Purpose: Compare ligand behavior across multiple viral protein structures.
- What the user can do:
  - Select a ligand or synonym.
  - View all mapped PDB-chain-residue contexts for that ligand.
  - Select multiple structures and generate comparison charts.
  - View the ligand’s 2D SVG.
  - Push the ligand into PROTAC Builder using the “Use ligand in PROTAC Builder” controls.
- Main controls/forms:
  - Ligand selector.
  - Multi-selection checkbox grid of PDB contexts.
  - “Generate comparison charts” button.
  - Ligand modal.
  - Plot carousel with:
    - Interaction type distribution
    - Interaction distance profile
    - Stacked interaction counts
    - Atom-level interaction burden
- Backend routes/APIs called:
  - `/get_ligands_with_synonyms`
  - `/get_pdb_mapping/<ligand_code>`
  - `/get_smiles_svg/<ligand_id>`
  - `/compare_ligand_interactions`
  - `/drugapp/` HEAD check for chatbot availability
- Database tables used:
  - `ligand_synonyms`
  - `Ligand_Atoms_Smiles`
  - `SMILES_MAP_PDB`
  - `Arpeggio_Contacts_Data`
- Output/export behavior:
  - Chart panels are rendered client-side with Plotly from JSON returned by `/compare_ligand_interactions`.
  - No direct CSV export is defined on this page, but the visual outputs support manual comparative review and downstream handoff.
- Scientific interpretation:
  - Good for determining whether a ligand’s interaction pattern is conserved or context-specific across structures.
  - Atom-level interaction burden can help prioritize atoms that are less interaction-critical and therefore more likely to tolerate modification.
- Known limitations:
  - Mapping quality depends on `SMILES_MAP_PDB`.
  - Color scaling and chart density become harder to parse when many structures are selected.
- Maintenance notes:
  - This page is a major source of warhead triage intuition even though it is not the formal PROTACability dashboard.

### PROTACability Assessment
- Route: `/protacability_page`
- Template: `templates/protacability_assessment.html`
- Purpose: Multi-view target triage dashboard for degrader-aware interpretation of viral protein structures.
- What the user can do:
  - Browse target-level, protein-level, structure-level, and chain-level evidence.
  - Filter on virus, protein type, scores, tiers, ligand context, and evidence flags.
  - Open detail modals with tabular context and NGL-based 3D viewing.
  - Export filtered results or full imported PROTACability tables as CSV.
- Main controls/forms:
  - Views:
    - Target Browser
    - Protein Summary
    - Structure Summary
    - Chain Details
  - Filters:
    - Virus
    - Protein type
    - Priority tier
    - Degrader readiness tier
    - Warhead linkability tier
    - PDB code
    - Minimum proxy score
    - Minimum degrader readiness score
    - Minimum warhead linkability score
    - Candidate ligand resname
    - Ligand presence
    - Ligand context
    - Evidence level
    - SMILES source
    - Collapse overlapping protein labels
    - Has exposed lysine
    - Has ligand-proximal exposed lysine
    - Has candidate linker atoms
    - Has solvent-exposed ligand atoms
    - Has mapped atoms
    - Has valid RDKit SMILES
    - Has exposed target lysines
  - Toolbar:
    - Sort selector
    - Page size selector
    - Apply filters
    - Reset
    - Export selector and CSV export button
- Backend routes/APIs called:
  - `/api/protacability/search`
  - `/api/protacability/filter_options`
  - `/api/protacability/filters`
  - `/api/protacability/detail/<pdb_code>/<chain_id>`
  - `/api/protacability/structure_detail/<pdb_code>`
  - `/api/protacability/protein_detail`
  - `/api/protacability/target_detail`
  - `/api/protacability/export`
  - `/api/coordinates/<pdb_code>.pdb`
  - `/api/ligand_instance_sdf/<pdb_code>/<ligand_code>.sdf`
  - `/api/ligand_instance_sdf_url/<pdb_code>/<ligand_code>`
  - Optional debug endpoints for ligand mapping and coordinate resolution
- Database tables used:
  - Core:
    - `protacability_assessment`
    - `protacability_lysine_proximity`
    - `protacability_ligand_inventory`
  - Optional enrichment:
    - `protacability_warhead_linkability`
    - `protacability_degrader_readiness`
  - Supporting structure delivery:
    - local coordinate files under `PDB_FILES/`
    - `static/coordinate_cache/`
    - `static/ligand_sdf_cache/`
- Output/export behavior:
  - Export filtered current-view results as `protacability_<view>_filtered.csv`.
  - Export full imported tables when a specific export table is chosen from the selector.
  - Open detail modal with NGL viewer and a “Use as Warhead” button for downstream design flow.
- Scientific interpretation:
  - The dashboard explicitly separates five evidence concepts:
    1. Overall Degrader Readiness
    2. Protein Structural Priority
    3. Warhead Linkability
    4. Target Lysine Accessibility
    5. Ternary Geometry Cue
  - Use the following framing when interpreting the outputs:
    - Warhead linkability evaluates whether the bound ligand contains solvent-exposed, chemically modifiable atoms that may tolerate linker attachment while preserving binding.
    - Target lysine accessibility is evaluated separately as a ubiquitination-readiness cue.
    - Ternary geometry cues are hypothesis-generating only and should not be interpreted as proof of productive degradation.
  - Avoid older shorthand such as “lysine-linker site,” “linker attaches to target lysine,” or “lysine linkability.”
- View descriptions:
  - Target Browser:
    - Best for target-level triage.
    - Collapses repeated structure contexts into virus/protein target group cards.
    - Surfaces overall degrader readiness while keeping warhead evidence separate from lysine accessibility and geometry cues.
  - Protein Summary:
    - Collapses repeated structures to make target-level comparison easier to scan.
    - Good for prioritizing protein classes rather than individual PDBs.
  - Structure Summary:
    - Focuses on representative PDB systems and summarized chain evidence.
    - Helpful for identifying concrete structure exemplars to carry into downstream visualization or chemistry.
  - Chain Details:
    - Exposes raw chain-level evidence and the richest feature context.
    - Best for checking mapped atoms, candidate linker atoms, and per-chain readiness fields.
- Data sources:
  - Existing imported tables:
    - `protacability_assessment`
    - `protacability_lysine_proximity`
    - `protacability_ligand_inventory`
  - New enrichment tables:
    - `protacability_warhead_linkability`
    - `protacability_degrader_readiness`
  - Generated CSVs:
    - `PDB_FILES/PROTACability_Assessment.csv`
    - `PDB_FILES/PROTACability_Lysine_Ligand_Proximity.csv`
    - `PDB_FILES/PROTACability_Ligand_Inventory.csv`
    - `PDB_FILES/PROTACability_Warhead_Linkability.csv`
    - `PDB_FILES/PROTACability_Degrader_Readiness.csv`
- Enrichment script:
  - `01_PROTACability_Warhead_Linkability_Enrichment.py`
  - Uses `viral_data.db` when available, plus optional `Components-smiles-stereo-oe.smi` fallback for component SMILES coverage.
  - Uses multiprocessing for faster ligand-centered scoring.
  - Does not prove experimental degradation.
  - Produces linkability/readiness evidence tables and a failure log.
- Known limitations:
  - Results are heuristic and data-dependent.
  - Candidate linker atom visualization is prepared in payloads and partly surfaced in detail views, but exact 3D highlighting remains an area for future refinement.
  - Chain duplication and overlapping labels are handled with collapse logic, but interpretation still requires judgment.
  - Glycan-only or common-buffer ligands are intentionally not treated as strong warhead evidence.
- Maintenance notes:
  - This is the most important page to keep consistent with manuscript language and downstream chemistry interpretation.
  - If enrichment tables are missing, the page still runs but with reduced evidence depth.

### PROTAC Builder
- Route:
  - External links from nav, footer, Home, Ligand Comparison, display pages, and legacy redirect `/copy`.
  - Configurable base URL from `PROTAC_BUILDER_EXTERNAL_URL`, defaulting to `https://protacbuilder.com`.
- Template:
  - Not hosted in this repository as a full internal module.
- Purpose:
  - Downstream modular degrader design workflow that receives promising ligands/warheads from V-LiSEMOD.
- What the user can do:
  - Open the standalone builder from multiple handoff points.
  - Use a ligand as a candidate warhead in downstream design.
- Main controls/forms:
  - In V-LiSEMOD, this appears as external-link buttons and “Use as Warhead” style controls.
- Backend routes/APIs called:
  - `/copy` and `/copy/<path:legacy_path>` redirect to the external builder URL.
  - Deprecated legacy routes such as `/run-drug-analysis`, `/download/<filename>`, `/process`, and `/ligase_details` now return deprecation JSON pointing users to PROTAC Builder.
- Database tables used:
  - None directly inside this repo’s route layer for the external builder itself.
- Output/export behavior:
  - Handoff only from the V-LiSEMOD perspective.
- Scientific interpretation:
  - V-LiSEMOD is the upstream triage and evidence workspace; PROTAC Builder is the downstream design environment.
- Known limitations:
  - TODO: Internal builder workflows such as recruiter/linker property filtering are not verifiable from this repository alone because the main builder implementation lives outside this codebase.
- Maintenance notes:
  - Keep link destinations and deprecation redirects aligned with the currently supported external builder deployment.

### Drug GPT / BioGPT
- Route:
  - `/drugapp/` and `/drugapp/query` when enabled.
  - Friendly disabled responses when not enabled.
- Template:
  - `templates/DRUGindex.html` when enabled.
- Purpose:
  - Local LLM-assisted molecule support workspace, separate from the core curated structure workflows.
- What the user can do:
  - In enabled mode, ask free-text questions and receive model responses with optional PubChem enrichment.
  - In default public/lab-safe mode, see that the module is disabled.
- Main controls/forms:
  - Chat-style prompt box in the Drug GPT page.
- Backend routes/APIs called:
  - Implemented in `DRUGapp.py` via Flask blueprint `dp`.
  - Loaded by `app.py` only when `ENABLE_DRUG_GPT=1`.
- Database tables used:
  - None required for core model inference.
  - Optional PubChem enrichment is external and API-based.
- Output/export behavior:
  - JSON chat responses from `/drugapp/query`.
- Scientific interpretation:
  - This is not a required part of the current V-LiSEMOD scientific workflow.
- Known limitations:
  - Hidden/disabled for current public/lab deployment.
  - Local LLM/GPU loading should be disabled by default.
  - Model downloads are expensive and should not happen automatically in lightweight deployments.
- Maintenance notes:
  - Code still exists and can be re-enabled later after domain-specific model work is ready.
  - Do not rely on Drug GPT for current project functionality.

## 6. Database Overview
The primary application database is `viral_data.db` (SQLite). The auth/session database is `users.db`.

### Core scientific database inventory
| Table name | Purpose | Key columns | Used by pages/modules | Notes |
|---|---|---|---|---|
| `Arpeggio_Contacts_Data` | Atom-level interaction/contact records between ligand atoms and residues | `virus_name`, `pdb_id`, `ligand`, `ligand_id`, `chain`, `Contact`, `Distance`, `exact_atom`, `atom_id`, `residue`, `residue_number` | Ligand Indexer, Ligand Comparison, export bundles, enrichment support | Core interaction evidence layer. |
| `Covalent_Noncovalent` | Annotates inhibitor type as covalent or noncovalent | `virus_name`, `pdb_id`, `ligand`, `ligand_id`, `chain`, `Inhibitor_Type` | Not obviously surfaced in current main pages | Present in DB; may support future classification or analysis. |
| `Functional_GROUPED` | Ligand SMILES and grouped functional group annotations | `virus_name`, `pdb_id`, `ligand`, `smiles`, `functional_groups` | Home ligand-image generation, exports | Used to generate 2D ligand display content. |
| `Functional_Group_Atoms` | Functional-group atom-level annotation | `virus_name`, `pdb_id`, `ligand`, `chain`, `functional_group`, `atom_id`, `exact_atom` | Home PyMOL session generation, exports | Enables functional group object generation. |
| `Ligand_Arp_Diagram` | Ligand occurrence mapping by PDB/chain/residue | `virus_name`, `pdb_id`, `ligand`, `chain`, `ligand_id` | Home, Ligand Indexer, exports, display pages | Primary ligand occurrence lookup table. |
| `Ligand_Atoms_Smiles` | Ligand mapping plus representative SMILES and MW | `virus_name`, `pdb_id`, `ligand`, `chain`, `ligand_id`, `smiles`, `molecular_weight` | Ligand Comparison, Protein Query exports, ligand lookup APIs | Main ligand structure lookup table. |
| `RUPLEY_SASA_DATA` | Shrake-Rupley solvent-accessible surface area results for ligand atoms | `virus_name`, `pdb_id`, `ligand`, `chain`, `exact_atom`, `atom_id`, `SASA_Area` | Home, display pages, export bundles, comparison support | Critical for solvent-exposed atom visualization. |
| `SMILES_MAP_PDB` | Atom mapping between PDB atoms and SMILES atom indices | `virus_name`, `pdb_id`, `ligand`, `chain`, `exact_atom`, `atom_id`, `smiles_atom_index` | Ligand Comparison, display pages, PROTACability enrichment, exports | Key bridge between structure-space and chemistry-space. |
| `Virus_Proteins` | Viral structure classification by virus and protein label | `virus_name`, `pdb_id`, `protein` | Protein Query | Central target-centric lookup table. |
| `distal_atoms` | Distal atom set for PyMOL annotation | `virus_name`, `pdb_id`, `ligand`, `chain`, `atom_id`, `exact_atom` | Home PyMOL session generation | Present in core builder, not emphasized in current UI copy. |
| `ligand_atoms` | Raw ligand atom coordinates and types | `virus_name`, `pdb_id`, `ligand`, `chain`, `atom_id`, `exact_atom`, `x`, `y`, `z` | Home, ligand lists, PyMOL generation, exports | Fundamental ligand coordinate table. |
| `ligand_synonyms` | Ligand-to-synonym normalization | `ligand`, `synonym` | Protein Query, Ligand Indexer, Ligand Comparison | App code uses mixed-case `Ligand_Synonyms` in places; SQLite resolves it case-insensitively. |
| `ligand_water_distances` | Ligand-to-water distance records | `virus_name`, `pdb_id`, `ligand`, `chain`, `atom_id`, `water_chain`, `distance` | Not obviously surfaced in current main pages | Useful candidate support table for hydration analysis. |
| `receptor_binding_pocket` | Binding pocket residue and atom records | `virus_name`, `pdb_id`, `residue`, `residue_chain`, `residue_number`, `residue_atom` | Home PyMOL session generation, exports | Used to define yellow binding-pocket surfaces. |
| `solvent_exposed_atoms` | Legacy/parallel solvent-exposed atom coordinate records | `virus_name`, `pdb_id`, `ligand`, `chain`, `atom_id`, `exact_atom` | Not primary in current UI; export/background support | Distinct from `RUPLEY_SASA_DATA`; verify role before refactoring. |

### PROTACability tables
| Table name | Purpose | Key columns | Used by pages/modules | Notes |
|---|---|---|---|---|
| `protacability_assessment` | Core structural-priority assessment output from `00_Protein_Expansion3.py` | `virus_name`, `protein_type`, `pdb_code`, `chain_id`, `candidate_ligand_resnames`, `exposed_lys_count`, `protacability_proxy_score`, `protacability_tier` | PROTACability dashboard, Protein Query exports | Base chain-level triage layer. |
| `protacability_lysine_proximity` | Lysine accessibility/proximity evidence | `virus_name`, `protein_type`, `pdb_code`, `chain_id`, `lys_residue_id`, `is_surface_exposed`, `nearest_ligand_resname`, `nearest_ligand_distance_a` | PROTACability detail views, exports | Supports target lysine accessibility interpretation. |
| `protacability_ligand_inventory` | Ligand inventory per structure/model | `virus_name`, `protein_type`, `pdb_code`, `ligand_resname`, `ligand_chain`, `ligand_residue_id`, `ligand_atom_count` | PROTACability detail views, exports, enrichment inputs | Ligand context inventory used for 3D selection. |
| `protacability_warhead_linkability` | Ligand-centered warhead evidence layer | `virus_name`, `protein_type`, `pdb_code`, `ligand_resname`, `ligand_chain`, `candidate_linker_atom_count`, `warhead_linkability_score`, `warhead_linkability_tier`, `smiles_source` | PROTACability dashboard, exports | Optional enrichment table imported after enrichment run. |
| `protacability_degrader_readiness` | Combined degrader-readiness summary layer | `virus_name`, `protein_type`, `pdb_code`, `chain_id`, `degrader_design_readiness_score`, `degrader_design_readiness_tier`, `evidence_level`, `best_linker_geometry_class` | PROTACability dashboard, exports | Optional enrichment table imported after enrichment run. |

### Auth/session database inventory
| Table name | Purpose | Key columns | Used by pages/modules | Notes |
|---|---|---|---|---|
| `users` in `users.db` | Basic user records | `id`, `email`, `password_hash`, `role`, `is_active` | Login scaffold | No secrets should be documented or exported. |
| `sessions` in `users.db` | Session concurrency tracking | `sid`, `user_id`, `created_at`, `last_seen`, `user_agent`, `ip` | Login scaffold, session management | Controlled by `MAX_SESSIONS_PER_USER` and `SESSION_CONCURRENCY_MODE`. |

## 7. File/Folder Structure
| Path | Role | Notes |
|---|---|---|
| `app.py` | Main Flask application | Contains routes, exports, auth scaffold, PROTACability APIs, coordinate/SDF serving, and deployment feature flags. |
| `DRUGapp.py` | Optional Drug GPT blueprint | Hidden/disabled by default; loads local LLM only when enabled. |
| `templates/` | Jinja templates | Main pages, modals, output pages, shared shell. |
| `templates/base.html` | Shared layout | Global nav, footer, dark/light shell, theme toggle. |
| `templates/_tutorial_flyover.html` | Shared tutorial/help banner | Reused on multiple pages. |
| `static/` | Static web assets | CSS, JS, images, generated charts, caches. |
| `static/js/` | Client-side behavior | Includes `scripts.js`, `theme-toggle.js`, `ngl_viewer_helpers.js`. |
| `static/css/` | Styling | Shared styles and theme variables. |
| `static/images/` | App image assets | Branding and tutorial/support imagery. |
| `static/ligand_images/` | Generated ligand SVG assets | Do not commit generated outputs. |
| `static/charts/` | Generated chart outputs | Ignored/generated. |
| `static/coordinate_cache/` | Cached coordinate conversions | Ignored/generated. |
| `static/ligand_sdf_cache/` | Cached ligand SDF files | Ignored/generated. |
| `PDB_FILES/` | Large structural dataset and generated PROTACability CSVs | Ignored/generated; not for Git. |
| `output_csvs/` | Intermediate/generated CSV inputs | Ignored/generated. |
| `pml_sessions/` | Generated PyMOL session outputs | Ignored/generated. |
| `docs/` | Human documentation | Good place for this file and future user guides. |
| `TOOLS/import_protacability_data.py` | CSV-to-SQLite import utility | Imports/updates PROTACability tables and indexes. |
| `00_Protein_Expansion.py` | Structure download/manifest pipeline | Builds `PDB_FILES/` structure set and manifest. |
| `00_Protein_Expansion2.py` | Lysine/pI/SASA expansion | Produces `PDB_FILES/Lysine_ISO.csv`. |
| `00_Protein_Expansion3.py` | Core PROTACability-style assessment generator | Produces assessment, lysine proximity, ligand inventory CSVs. |
| `01_PROTACability_Warhead_Linkability_Enrichment.py` | Ligand-centered enrichment layer | Produces warhead linkability and degrader readiness CSVs. |
| `README.md` | Project-level readme | Includes deployment-default notes and environment toggles. |
| `requirements.txt` | Python dependency list | Includes Flask, pandas, RDKit, matplotlib, torch, transformers, biopython, etc. |

### Folders/files that should not be committed
These are either already ignored in `.gitignore` or should be treated as deployment-specific/generated assets:

- `PDB_FILES/`
- Generated PROTACability CSVs
- `viral_data.db`
- `users.db`
- `Components-smiles-stereo-oe.smi`
- `static/coordinate_cache/`
- `static/ligand_sdf_cache/`
- `static/charts/`
- `static/ligand_images/` when generated locally
- `output_csvs/`
- `pml_sessions/`
- `.env`
- model files and checkpoints

## 8. Deployment / Environment Notes
Default intended behavior for the current V-LiSEMOD deployment:

- The app should run without GPU by default.
- Drug GPT/local LLM loading should be disabled unless explicitly enabled.
- Large data files are ignored by Git and must be provisioned or regenerated separately.
- The scientific database (`viral_data.db`) must be provisioned separately from source code.
- PDB files and caches are regenerated or stored outside Git.

Implemented environment variables and relevant config:

| Variable | Purpose | Current status |
|---|---|---|
| `ENABLE_DRUG_GPT` | Enables Drug GPT blueprint/routes | Implemented in `app.py` |
| `ENABLE_LOCAL_LLM` | Allows local LLM configuration/model loading | Implemented in `app.py` and `DRUGapp.py` |
| `SHOW_DRUG_GPT_NAV` | Shows Drug GPT nav item | Implemented in `app.py` / `base.html` |
| `FLASK_SECRET_KEY` | Flask session secret | Implemented |
| `MAX_SESSIONS_PER_USER` | Max concurrent sessions per user | Implemented |
| `SESSION_CONCURRENCY_MODE` | Session handling policy such as `DENY` or `EVICT` | Implemented |
| `PROTAC_BUILDER_EXTERNAL_URL` | External builder base URL | Implemented |
| `MODEL_ID` / `LLM_MODEL_ID` | Override local model ID | Implemented for Drug GPT |
| `HF_TOKEN` / Hugging Face token variants | Optional authenticated model download | Implemented for Drug GPT |

Auth note:

- Login/session scaffolding exists and `users.db` stores users and sessions.
- However, the current `EXEMPT_PATH_PREFIXES` tuple in `app.py` includes `/`, which effectively exempts all routes from login enforcement because every path starts with `/`.
- For many lab/public deployments this may be intentional, but it should be treated as a current behavior to verify rather than an assumed security model.

Run command examples visible in the repo:

- `waitress-serve --listen=127.0.0.1:5002 app:app`
- `flask run`

## 9. Maintenance Checklist
### Updating and regenerating data
1. Update or regenerate source structure lists feeding `output_csvs/` and manifest creation.
2. Rebuild/download PDB/mmCIF content:
   - `python 00_Protein_Expansion.py`
3. Recompute lysine/pI/SASA summary outputs if needed:
   - `python 00_Protein_Expansion2.py`
4. Rebuild PROTACability-style core structural assessment tables:
   - `python 00_Protein_Expansion3.py --fresh --workers 8`
5. Recompute ligand-centered warhead/readiness enrichment:
   - `python 01_PROTACability_Warhead_Linkability_Enrichment.py --workers 8 --component-smiles Components-smiles-stereo-oe.smi`
6. Import generated PROTACability CSVs into SQLite:
   - `python TOOLS/import_protacability_data.py`

### Verifying the app
1. Start locally:
   - `waitress-serve --listen=127.0.0.1:5002 app:app`
2. Check route health:
   - Visit `/`
   - Visit `/about`
   - Visit `/query_protein_virus_page`
   - Visit `/ligand_indexer`
   - Visit `/compare_ligands`
   - Visit `/protacability_page`
   - Visit `/healthz`
3. Confirm that:
   - Navigation renders correctly.
   - Protein Query exports still download.
   - Ligand Comparison charts still generate.
   - PROTACability filters populate and detail modal opens.
   - Drug GPT remains disabled unless intentionally enabled.

### Avoiding bad commits
1. Review `.gitignore` before staging.
2. Check status safely:
   - `git status`
3. Avoid adding:
   - databases
   - `PDB_FILES/`
   - caches
   - generated charts
   - generated ligand images
   - model artifacts
   - `.env`

### Notes on specific scripts
- `00_Protein_Expansion3.py` supports:
  - `--fresh`
  - `--resume`
  - `--skip-failures`
  - `--workers`
  - `--limit`
  - `--sasa-n-points`
  - `--probe-radius`
- `01_PROTACability_Warhead_Linkability_Enrichment.py` supports:
  - `--workers`
  - `--serial`
  - `--limit`
  - `--progress-every`
  - `--component-smiles`
  - path overrides for inputs and outputs

## 10. Known Limitations
- PROTACability and degrader readiness are heuristic, not experimentally validated.
- Ligand linkability depends on SMILES, atom mapping, SASA, and contact completeness.
- Component SMILES fallback helps coverage but may be lower confidence than curated ligand-specific data.
- Glycans and common crystallization/buffer ligands are intentionally not treated as strong warhead evidence.
- Some PDB structures may have chain duplication or overlapping protein labels.
- Some ligand mappings may be incomplete or require alias handling.
- NGL highlighting infrastructure exists, but future work is still needed for polished candidate-linker-atom visualization in 3D.
- Drug GPT is disabled for now and should not be relied on in current workflows.
- Auth/session scaffolding exists, but current route exemption behavior effectively leaves the app public unless changed.

## 11. Future Work
- Improve candidate linker atom visualization in 3D detail views.
- Add better manuscript-quality exports and reporting.
- Add a model-trained AI assistant only after a domain-specific model is ready and deployment-safe.
- Integrate PROTAC Builder outputs more tightly with V-LiSEMOD context and return paths.
- Add automated tests for routes, APIs, database expectations, and import scripts.
- Add a validation report after each large regeneration/import run.
- Expand user-facing tutorial and help documentation beyond the current tutorial banner pattern.

## 12. Quick Start for a New Lab User
1. Open the app in your browser and start at Home.
2. If you want a target-centric workflow, go to Protein Query.
3. Select a virus and protein type to narrow to relevant viral targets.
4. If you want a ligand-first workflow, use Ligand Indexer or Ligand Comparison.
5. For degrader-oriented triage, open PROTACability Assessment and start in Target Browser.
6. Use the filtered views and evidence tiers to identify promising targets or ligands.
7. Export data when you need a structure bundle, CSVs, or a PyMOL review package.
8. Do not overinterpret a high readiness score as proof of degradability; treat it as a structured starting point for discussion and follow-up.

## 13. Quick Start for a Maintainer
1. Install dependencies from `requirements.txt` or the project conda environment.
2. Provision `viral_data.db` outside Git.
3. Start locally with:
   - `waitress-serve --listen=127.0.0.1:5002 app:app`
4. Verify the major pages and export routes manually.
5. When source data changes, rerun:
   - `python 00_Protein_Expansion.py`
   - `python 00_Protein_Expansion2.py`
   - `python 00_Protein_Expansion3.py --fresh --workers 8`
   - `python 01_PROTACability_Warhead_Linkability_Enrichment.py --workers 8 --component-smiles Components-smiles-stereo-oe.smi`
   - `python TOOLS/import_protacability_data.py`
6. Avoid committing:
   - `PDB_FILES/`
   - `viral_data.db`
   - `users.db`
   - caches
   - generated exports
   - `.env`
7. To re-enable Drug GPT later:
   - set `ENABLE_DRUG_GPT=1`
   - set `ENABLE_LOCAL_LLM=1` only if model loading is truly desired
   - optionally set `SHOW_DRUG_GPT_NAV=1`
   - verify the Home template and any chatbot stubs still behave correctly under the chosen deployment mode
