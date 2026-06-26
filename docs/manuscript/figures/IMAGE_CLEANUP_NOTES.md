# Image Cleanup Notes

## Images needing crop

- `raw/fig01_vlisemod_home_full.png`
- `raw/fig03_warhead_browse_full.png`
- `raw/fig04_vlisemod_home_full.png`
- `raw/fig04_vlisemod_protacability_full.png`
- `raw/fig05_protacbuilder_builder_full.png`
- `raw/fig06_e3ligandalyzer_home_full.png`
- `raw/fig06_e3ligandalyzer_explorer_full.png`
- `raw/fig01_pymacs_github_full.png`
- `raw/fig06_pymacs_github_full.png`
- `raw/fig07_warhead_api_docs_full.png`

## Images needing annotation

- `annotated/fig01_ecosystem_overview_draft.png`
  - add panel labels and possibly arrows or short tool-role labels
- `cropped/fig04_vlisemod_protacability_filters_results.png`
  - likely benefits from subtle label callouts for heuristic triage fields
- `cropped/fig03_warhead_exposure_result_panel.png`
  - likely benefits from callouts clarifying that the panel is a public browse view, not experimental validation

## Images needing redaction

- No obvious credentials, cookies, or local file paths were seen in the spot-checked captures.
- Recheck GitHub screenshots before publication to confirm no personal sign-in state or UI profile data appears in the final exported panel.

## Screenshots with too much whitespace or repeated content

- `raw/fig01_vlisemod_home_full.png`
- `raw/fig03_warhead_browse_full.png`
- `raw/fig04_vlisemod_home_full.png`
- `raw/fig05_protacbuilder_builder_full.png`
- `raw/fig06_e3ligandalyzer_home_full.png`
- `raw/fig06_e3ligandalyzer_explorer_full.png`

These pages repeated content in full-page mode because of sticky or dynamic layout behavior.

## Screenshots with tiny text

- `raw/fig07_warhead_api_docs_full.png`
- `raw/fig01_pymacs_github_full.png`
- `raw/fig06_pymacs_github_full.png`
- `raw/fig03_warhead_browse_full.png`

## Pages that should be recaptured after UI cleanup or local app access

- AutoDock-Vina PrepServer figure slot once the correct local repository or running app is available
- E3 Ligandalyzer pages if a cleaner non-overlay capture is desired
- PROTAC Builder builder workspace if a post-modal public state can be captured safely without interacting deeply

## Suggested final composite layout

- Figure 1: five-panel ecosystem montage using the current draft composite as a storyboard base
- Figure 3: three-panel Warhead Hunter workflow composite using hunter, examples, and results-library crops
- Figure 4: five-panel V-LiSEMOD composite using structure explorer, Protein Query, Ligand Indexer, Ligand Comparison, and PROTACability crops
- Figure 5: two-panel PROTAC Builder continuation composite
- Figure 6: two-panel E3 Ligandalyzer plus PyMACS companion-context composite
- Figure 7: one screenshot plus a small inset summarizing API evidence files
