# Manuscript Draft Notes

## Source Files Used

- `README.md`
- `docs/APP_GUIDE.md`
- `docs/DATABASE.md`
- `docs/DEVELOPER_NOTES.md`
- `docs/PROTACABILITY.md`
- `docs/MANUSCRIPT_OUTLINE.md`
- `docs/manuscript/MANUSCRIPT_PLAN.md`
- `docs/manuscript/CLAIMS_AND_LIMITATIONS_MATRIX.md`
- `docs/manuscript/FIGURE_AND_TABLE_PLAN.md`
- `docs/manuscript/VALIDATION_AND_REPRODUCIBILITY_PLAN.md`
- `docs/manuscript/figures/FIGURE_INDEX.md`
- `docs/manuscript/figures/FIGURE_CAPTIONS.md`
- `docs/manuscript/figures/FIGURE_STORYBOARD.md`
- `docs/manuscript/figures/SCREENSHOT_CAPTURE_NOTES.md`
- `docs/manuscript/figures/IMAGE_CLEANUP_NOTES.md`
- `docs/manuscript/figures/TOOL_URL_CHECKS.md`
- `docs/manuscript/figures/api_evidence/README.md`

## Figures Inserted

- Figure 1: docs/manuscript/figures/annotated/fig01_ecosystem_overview_draft.png
- Figure 3: docs/manuscript/figures/cropped/fig03_warhead_hunter_panel.png
- Figure 3: docs/manuscript/figures/cropped/fig03_warhead_exposure_result_panel.png
- Figure 4: docs/manuscript/figures/cropped/fig04_vlisemod_home_panel.png
- Figure 4: docs/manuscript/figures/cropped/fig04_vlisemod_protein_query_panel.png
- Figure 4: docs/manuscript/figures/cropped/fig04_vlisemod_ligand_indexer_panel.png
- Figure 4: docs/manuscript/figures/cropped/fig04_vlisemod_compare_ligands_panel.png
- Figure 4: docs/manuscript/figures/cropped/fig04_vlisemod_protacability_filters_results.png
- Figure 5: docs/manuscript/figures/cropped/fig05_protacbuilder_home_panel.png
- Figure 5: docs/manuscript/figures/cropped/fig05_protacbuilder_design_workspace.png
- Figure 6: docs/manuscript/figures/cropped/fig06_e3_pymacs_companion_panels.png
- Figure 7: docs/manuscript/figures/raw/fig07_warhead_api_docs_full.png

## Figures Omitted or Missing

- Figure 2: Image pending: regenerate this panel once the correct local AutoDock-Vina PrepServer application or repository is available.

## Tables Included

- Platform capabilities
- Data-layer summary
- PROTACability interpretation guardrails
- Module-to-user-question map
- Limitations and future work
- Reproducibility checklist

## Claims Intentionally Avoided

- Experimentally validated degradation prediction
- Guaranteed PROTAC design success
- Productive ternary-complex prediction
- Automated medicinal chemistry decision-making
- Exhaustive viral ligand or target coverage
- Full technical unification across companion tools
- Production-scale API guarantees or service-level commitments

## Citation Placeholders Needing Follow-up

- [CITATION NEEDED: RCSB/PDB or relevant structural database]
- [CITATION NEEDED: Flask/software framework]
- [CITATION NEEDED: SQLite or database layer]
- [CITATION NEEDED: PyMOL if export workflow is discussed]
- [CITATION NEEDED: Arpeggio/contact-analysis method if cited]
- [CITATION NEEDED: solvent-accessible surface area method]
- [CITATION NEEDED: PROTAC review]
- [CITATION NEEDED: targeted protein degradation review]
- [CITATION NEEDED: viral structural bioinformatics or antiviral ligand-design context]

## Validation Notes

- The DOCX should be validated after regeneration with `unzip -t docs/manuscript/drafts/V-LiSEMOD_manuscript_draft.docx`.
- The DOCX should be loaded with `python-docx` to confirm paragraphs, tables, and inline images are present.
- A render attempt with the bundled `render_docx.py` may fail on this workstation if headless LibreOffice cannot load `/opt/homebrew/opt/little-cms2/lib/liblcms2.2.dylib`; that is a local renderer dependency blocker, not evidence that the DOCX package is corrupt.
- If LibreOffice is repaired or available elsewhere, rerun DOCX-to-PNG rendering and visually inspect all pages before journal submission.

## Manual Review Checklist

- Add real author names, affiliations, contact details, funding, acknowledgments, conflicts, and ethics statements as appropriate.
- Confirm final repository URL, public URLs, license, and data availability language.
- Replace citation placeholders with verified references.
- Recheck all figure panels for publication readiness, overlays, modals, tiny text, and private UI state.
- Regenerate the pending AutoDock-Vina PrepServer panel from the correct application context.
- Re-run local app validation and endpoint checks against the final manuscript version.
- Confirm PROTACability wording remains explicitly heuristic throughout the manuscript.
- Decide whether the article target is a software article, webserver paper, application note, methods article, or resource article.
