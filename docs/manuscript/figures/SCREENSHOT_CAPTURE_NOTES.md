# Screenshot Capture Notes

## Capture setup

- Capture date: 2026-06-11
- Viewport: `1440x1000`
- Source type: public live websites plus local repository inspection
- Browser state: no credentials, no private session content, no local file paths intentionally displayed

## URLs visited

- `https://warheadhunter.com/`
- `https://warheadhunter.com/hunter`
- `https://warheadhunter.com/browse`
- `https://warheadhunter.com/examples`
- `https://warheadhunter.com/api-docs`
- `https://vlisemod.com/`
- `https://vlisemod.com/about`
- `https://vlisemod.com/query_protein_virus_page`
- `https://vlisemod.com/ligand_indexer`
- `https://vlisemod.com/compare_ligands`
- `https://vlisemod.com/protacability_page`
- `https://vlisemod.com/healthz`
- `https://vlisemod.com/drugapp/`
- `https://protacbuilder.com/`
- `https://protacbuilder.com/builder`
- `https://e3ligandalyzer.com/`
- `https://e3ligandalyzer.com/explorer`
- `https://github.com/schurerlab/Pymacs`

## Screenshots captured

### Raw

- `raw/fig01_warhead_home_full.png`
- `raw/fig01_vlisemod_home_full.png`
- `raw/fig01_protacbuilder_home_full.png`
- `raw/fig01_e3ligandalyzer_home_full.png`
- `raw/fig01_pymacs_github_full.png`
- `raw/fig03_warhead_hunter_full.png`
- `raw/fig03_warhead_examples_full.png`
- `raw/fig03_warhead_browse_full.png`
- `raw/fig04_vlisemod_home_full.png`
- `raw/fig04_vlisemod_about_full.png`
- `raw/fig04_vlisemod_protein_query_full.png`
- `raw/fig04_vlisemod_ligand_indexer_full.png`
- `raw/fig04_vlisemod_compare_ligands_full.png`
- `raw/fig04_vlisemod_protacability_full.png`
- `raw/fig05_protacbuilder_home_full.png`
- `raw/fig05_protacbuilder_builder_full.png`
- `raw/fig06_e3ligandalyzer_home_full.png`
- `raw/fig06_e3ligandalyzer_explorer_full.png`
- `raw/fig06_pymacs_github_full.png`
- `raw/fig07_warhead_api_docs_full.png`

### Cropped

- `cropped/fig03_warhead_exposure_result_panel.png`
- `cropped/fig03_warhead_hunter_panel.png`
- `cropped/fig04_vlisemod_home_panel.png`
- `cropped/fig04_vlisemod_protein_query_panel.png`
- `cropped/fig04_vlisemod_ligand_indexer_panel.png`
- `cropped/fig04_vlisemod_compare_ligands_panel.png`
- `cropped/fig04_vlisemod_protacability_filters_results.png`
- `cropped/fig05_protacbuilder_home_panel.png`
- `cropped/fig05_protacbuilder_design_workspace.png`
- `cropped/fig06_e3ligandalyzer_module_panel.png`
- `cropped/fig06_e3ligandalyzer_explorer_panel.png`
- `cropped/fig06_pymacs_repo_panel.png`
- `cropped/fig06_e3_pymacs_companion_panels.png`

### Annotated or composite draft

- `annotated/fig01_ecosystem_overview_draft.png`

## API evidence captured

- `api_evidence/warhead_api_health.json`
- `api_evidence/warhead_api_manifest.json`
- `api_evidence/warhead_api_examples.json`
- `api_evidence/vlisemod_healthz.txt`

## Failed or unavailable pages

- `https://vlisemod.com/drugapp/`
  - returned `503`
  - documented as optional-disabled behavior rather than a missing public workflow

- Local AutoDock-Vina PrepServer routes suggested by the prompt
  - not captured
  - current workspace does not match that application

- `https://e3ligandalyzer.com/modules`
  - returned `404`

- `https://protacbuilder.com/resources`
  - returned `404`

- `https://protacbuilder.com/science`
  - returned `404`

## Rendering issues observed

- Several long full-page screenshots repeated page sections because sticky or dynamic layouts were duplicated in full-page capture mode.
- E3 Ligandalyzer screenshots captured with a visible modules dropdown overlay.
- PROTAC Builder builder screenshot opened with a welcome modal, which is still useful as a public workflow snapshot but may need recapture if a cleaner builder workspace panel is preferred.
- GitHub full-page screenshots are extremely tall and should be used only through cropped panels.

## Skipped items

- No credentialed routes were visited.
- No private admin panels were visited.
- No heavy jobs were submitted.
- No local structures or ligand libraries were uploaded to public tools.
