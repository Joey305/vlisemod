# Tool URL Checks

Capture session date: 2026-06-11  
Viewport target: `1440x1000`  
Capture mode: public live websites plus local repository inspection

## AutoDock-Vina PrepServer

| Tool | Base URL | Routes or pages checked | Observed behavior | Screenshot files | Notes |
|---|---|---|---|---|---|
| AutoDock-Vina PrepServer | local app expected by prompt | local repo pages, route hints from prompt | not captured from this workspace | none | This workspace is the V-LiSEMOD repository, not an AutoDock-Vina PrepServer repository. No matching local pages such as `/workflow`, `/modules`, or `/documentation` were found here, so the AutoDock figure slot remains pending a correct local app or separate repo. |

## Warhead Hunter

| Tool | Base URL | Routes or pages checked | Observed behavior | Screenshot files | Notes |
|---|---|---|---|---|---|
| Warhead Hunter | `https://warheadhunter.com` | `/`, `/hunter`, `/browse`, `/examples`, `/api-docs` | all returned `200` and rendered publicly | `raw/fig01_warhead_home_full.png`, `raw/fig03_warhead_hunter_full.png`, `raw/fig03_warhead_browse_full.png`, `raw/fig03_warhead_examples_full.png`, `raw/fig07_warhead_api_docs_full.png` | Good public coverage for ecosystem, workflow, results-library, examples, and API-docs panels. |

## V-LiSEMOD

| Tool | Base URL | Routes or pages checked | Observed behavior | Screenshot files | Notes |
|---|---|---|---|---|---|
| V-LiSEMOD | `https://vlisemod.com` | `/`, `/about`, `/query_protein_virus_page`, `/ligand_indexer`, `/compare_ligands`, `/protacability_page`, `/healthz`, `/drugapp/` | main public pages returned `200`; `/healthz` returned `200`; `/drugapp/` returned `503` | `raw/fig01_vlisemod_home_full.png`, `raw/fig04_vlisemod_home_full.png`, `raw/fig04_vlisemod_about_full.png`, `raw/fig04_vlisemod_protein_query_full.png`, `raw/fig04_vlisemod_ligand_indexer_full.png`, `raw/fig04_vlisemod_compare_ligands_full.png`, `raw/fig04_vlisemod_protacability_full.png` | `503` on `/drugapp/` is consistent with optional-disabled assistant behavior and was not used as a central figure. |

## PROTAC Builder

| Tool | Base URL | Routes or pages checked | Observed behavior | Screenshot files | Notes |
|---|---|---|---|---|---|
| PROTAC Builder | `https://protacbuilder.com` | `/`, `/builder`, plus checks for `/resources` and `/science` | `/` and `/builder` returned `200`; `/resources` and `/science` returned `404` | `raw/fig01_protacbuilder_home_full.png`, `raw/fig05_protacbuilder_home_full.png`, `raw/fig05_protacbuilder_builder_full.png` | Builder page is public and screenshot-safe, but only public no-input views were captured. No jobs were submitted. |

## E3 Ligandalyzer

| Tool | Base URL | Routes or pages checked | Observed behavior | Screenshot files | Notes |
|---|---|---|---|---|---|
| E3 Ligandalyzer | `https://e3ligandalyzer.com` | `/`, `/explorer`, `/ligases`, `/modules` | `/`, `/explorer`, and `/ligases` returned `200`; `/modules` returned `404` | `raw/fig01_e3ligandalyzer_home_full.png`, `raw/fig06_e3ligandalyzer_home_full.png`, `raw/fig06_e3ligandalyzer_explorer_full.png` | Homepage and explorer page are suitable for companion-context screenshots. |

## PyMACS

| Tool | Base URL | Routes or pages checked | Observed behavior | Screenshot files | Notes |
|---|---|---|---|---|---|
| PyMACS | `https://github.com/schurerlab/Pymacs` | repository landing page and `#readme` anchor check | returned `200` | `raw/fig01_pymacs_github_full.png`, `raw/fig06_pymacs_github_full.png` | Treated as a companion GitHub repository rather than a live web service. |
