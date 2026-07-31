# PROTACability Guide

## Core Interpretation

In V-LiSEMOD, PROTACability means a transparent, structure-guided prioritization framework for degrader-oriented follow-up. It combines ligand-centered and target-centered evidence into a practical triage layer for viral protein structures.

## What It Does Mean

- a way to separate promising ligand contexts from weaker ones
- a way to review solvent-exposed ligand atoms and potential linker-attachment cues
- a way to review exact ligand-instance attachment-site regions when attachment enrichment is available
- a way to evaluate whether target-side exposed lysines exist in the same structural context
- a way to rank structures for deeper design review

## What It Does Not Mean

- experimentally validated targeted degradation
- proof of productive ternary complex formation
- guaranteed warhead tractability or medicinal chemistry success
- a substitute for biochemical, cellular, or degradation assays

## Strong Interpretation Warning

PROTACability outputs are transparent structural-priority and design-readiness heuristics. They support hypothesis generation and triage, but they are not experimentally validated degradation predictions.

## Evidence Layers

### Warhead Linkability

- Focus: ligand-centered attachment opportunity
- Interprets: solvent-exposed mapped atoms, contact preservation, functional-group context, and chemically interpretable linkerability cues
- Best phrasing: "the bound ligand shows stronger or weaker structure-supported linker-attachment potential"

### Candidate Attachment Sites

- Focus: exact ligand-instance attachment regions generated from `attachment_v1_1`
- Interprets: candidate atom scores, region-level exposure, interaction cautions, and PDB atom serials for optional 3D highlighting
- Best phrasing: "this ligand instance has structure-supported candidate attachment-site regions"
- Boundary: these regions do not claim retained affinity, synthetic feasibility, ternary-complex formation, or degradation

### Target Lysine Accessibility

- Focus: target-side accessibility cues
- Interprets: presence and proximity of exposed lysines in the structural context
- Best phrasing: "the structure contains exposed lysine cues that may support degrader-oriented follow-up"

### Protein Structural Priority

- Focus: overall structure usefulness for follow-up
- Interprets: target context, ligand context, accessible lysine context, and aggregated evidence layers
- Best phrasing: "this structure ranks higher as a design-review starting point"

### Ternary Geometry Cue

- Focus: geometry-informed heuristic cue
- Interprets: ligand and lysine spatial context as a transparent approximation of whether the structure may be more compatible with future degrader exploration
- Best phrasing: "a hypothesis-generating geometry cue"

### Overall Degrader Readiness

- Focus: combined triage score
- Interprets: integrated heuristic readiness rather than proof of degradability
- Best phrasing: "overall degrader-readiness triage" or "combined structural-priority score"

## Preferred Interpretation Language

Use language such as:

- structure-supported linker-attachment opportunity
- transparent degrader-readiness heuristic
- ligand-centered warhead evidence
- target-side lysine accessibility cue
- prioritization for follow-up design review

## Wording To Avoid

Avoid language such as:

- predicted degrader
- validated degradation score
- guaranteed ternary complex
- proven degradability
- experimentally confirmed PROTACability

## Known Limitations

- limited by the quality and completeness of available co-crystal structures
- dependent on ligand mapping, SMILES coverage, and atom-level annotation quality
- attachment-site enrichment is only available for ligand instances that passed the graph/mapping eligibility checks
- lysine exposure is a structural cue, not a ubiquitination outcome
- geometry cues are simplified and do not model full ternary-complex behavior
- absence of a strong score does not prove a target is not degradable
- presence of a strong score does not guarantee successful degrader design

## Regeneration / Import Workflow

The app expects PROTACability CSV outputs in `PDB_FILES/` and imports them into `viral_data.db` using `TOOLS/import_protacability_data.py`.

Expected inputs:

- `PDB_FILES/PROTACability_Assessment.csv`
- `PDB_FILES/PROTACability_Lysine_Ligand_Proximity.csv`
- `PDB_FILES/PROTACability_Ligand_Inventory.csv`
- optional `PDB_FILES/PROTACability_Warhead_Linkability.csv`
- optional `PDB_FILES/PROTACability_Degrader_Readiness.csv`
- optional SQLite attachment enrichment tables: `protacability_attachment_analysis`, `protacability_attachment_regions`, and `protacability_attachment_atoms`

Operational note:

- the importer recreates table contents from the CSVs and adds indices for the main search surfaces

## Validation Checklist

Before presenting or exporting PROTACability results:

1. Confirm the relevant PROTACability tables are present in `viral_data.db`.
2. Confirm the dashboard no longer shows the "data not imported yet" message.
3. Spot-check a few structures for ligand presence, chain mapping, and exposed lysine context.
4. Spot-check at least one attachment-enriched ligand instance in the Candidate Attachment Sites detail section.
5. Verify that warhead/linkability and attachment-site language remains heuristic and non-claiming.
6. Verify that exported CSVs reflect the currently selected view or named legacy dataset.

## Export Behavior

The dashboard exposes CSV export for:

- the filtered current view
- legacy assessment-style exports
- lysine proximity
- ligand inventory
- warhead linkability
- degrader readiness
- attachment analysis
- attachment atoms
- attachment regions

Exports should be treated as analysis outputs for review, not as standalone evidence of biological activity.
