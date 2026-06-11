# Validation And Reproducibility Plan

## 1. Local Run Validation

1. Install dependencies from `requirements.txt` or `environment.yml`.
2. Start the application with `python app.py`.
3. Verify `/healthz`.
4. Verify the main pages:
   - `/`
   - `/query_protein_virus_page`
   - `/ligand_indexer`
   - `/compare_ligands`
   - `/protacability_page`
5. If the optional assistant is relevant, verify both disabled and enabled behavior for `/drugapp/`.

## 2. Data Validation

1. Verify `viral_data.db` exists when using local mode.
2. Confirm the expected core tables are present.
3. Confirm PROTACability tables are present:
   - `protacability_assessment`
   - `protacability_lysine_proximity`
   - `protacability_ligand_inventory`
   - `protacability_warhead_linkability`
   - `protacability_degrader_readiness`
4. Confirm generated cache and output folders are writable.
5. If RANDY mode will be described, verify the configured RANDY routes respond successfully.

## 3. Workflow Validation

### Structure Explorer

- select a virus, PDB, and ligand
- generate ligand imagery
- test PyMOL-oriented session generation
- confirm exposed-atom and functional-group overlays render when data exist

### Protein Query

- filter by virus and protein type
- verify PDB selection and export-linked dataset assembly
- verify PROTACability handoff when tables are present

### Ligand Indexer

- search a ligand by code or synonym
- verify mapped PDB and chain options
- generate an interaction chart

### Ligand Comparison

- compare one ligand across multiple structures
- verify chart generation and mapped-context behavior
- verify companion-tool handoff behavior

### PROTACability Dashboard

- verify filter loading
- verify search results across target, protein, structure, and chain views
- open detail views
- export filtered outputs

### PROTAC Builder Handoff

- verify the external handoff URL is correct for the deployment
- capture at least one representative handoff example

### PyMOL or Session Export

- confirm generated session files are created successfully
- verify outputs are written to expected runtime folders

## 4. Manuscript Evidence Checklist

- representative screenshots for each major app module
- at least one example target and ligand case
- example exported CSV outputs
- interaction chart examples
- ligand image examples
- PROTACability table and detail-view examples
- companion-tool handoff screenshots or URLs
- route and feature inventory for software description

## 5. Reproducibility Checklist

- exact app version or commit recorded
- environment variables documented
- database provisioning approach described
- data not committed publicly unless deliberately curated and redacted
- example workflow captured step by step
- limitations stated explicitly in the manuscript

## 6. Limitations Checklist

- co-crystal structure dependence
- database coverage dependence
- heuristic PROTACability interpretation
- no experimental degradation claims
- no guaranteed ternary-complex prediction
- optional local assistant is deployment-dependent
- generated outputs depend on local writable storage
- RANDY-backed deployment is environment-specific
