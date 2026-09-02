# V-LiSEMOD final release validation

* Validator: `final-release-validation-cif-v2.4`
* Database: `<author-build-directory>/viral_data_cif_v2.db`
* Final release validation: **PASS**
* Failed checks: 0
* Warnings: 0

This is a read-only release gate. It validates the frozen evidence generations and does not recalculate structural analyses.

## Schema

* **PASS** - table exists: structures: observed=1; expected=1
* **PASS** - table exists: ligands: observed=1; expected=1
* **PASS** - table exists: ligand_instances: observed=1; expected=1
* **PASS** - table exists: ligand_instance_atoms: observed=1; expected=1
* **PASS** - table exists: analysis_runs: observed=1; expected=1
* **PASS** - table exists: pipeline_failures: observed=1; expected=1
* **PASS** - table exists: mapping_remediation_queue: observed=1; expected=1
* **PASS** - table exists: ligand_mapping_runs: observed=1; expected=1
* **PASS** - table exists: ligand_smiles_atom_mapping: observed=1; expected=1
* **PASS** - table exists: ligand_sasa_atoms: observed=1; expected=1
* **PASS** - table exists: ligand_arpeggio_runs: observed=1; expected=1
* **PASS** - table exists: arpeggio_raw_contact_labels: observed=1; expected=1
* **PASS** - table exists: arpeggio_unique_atom_pairs: observed=1; expected=1
* **PASS** - table exists: protacability_ligand_inventory: observed=1; expected=1
* **PASS** - table exists: ligand_atom_geometry: observed=1; expected=1
* **PASS** - table exists: target_chain_geometry: observed=1; expected=1
* **PASS** - table exists: target_surface_lysines: observed=1; expected=1
* **PASS** - table exists: ligand_functional_group_matches: observed=1; expected=1
* **PASS** - table exists: ligand_functional_group_atoms: observed=1; expected=1
* **PASS** - table exists: ligand_functional_group_summary: observed=1; expected=1
* **PASS** - table exists: protacability_warhead_linkability: observed=1; expected=1
* **PASS** - table exists: protacability_target_context: observed=1; expected=1
* **PASS** - table exists: protacability_assessment: observed=1; expected=1
* **PASS** - table exists: protacability_degrader_readiness: observed=1; expected=1
* **PASS** - table exists: protacability_attachment_sites: observed=1; expected=1
* **PASS** - table exists: protacability_attachment_site_summary: observed=1; expected=1

## Compatibility

* **PASS** - view exists: v2_structure_context: observed=1; expected=1
* **PASS** - view exists: v2_ligand_context: observed=1; expected=1
* **PASS** - view exists: v2_ligand_atom_evidence: observed=1; expected=1
* **PASS** - view exists: v2_ligand_comparison_atom_contacts: observed=1; expected=1
* **PASS** - view exists: v2_all_chain_lysine_geometry: observed=1; expected=1
* **PASS** - view exists: v2_target_lysine_accessibility: observed=1; expected=1
* **PASS** - view exists: v2_protacability_target_context: observed=1; expected=1
* **PASS** - view exists: v2_protacability_best: observed=1; expected=1
* **PASS** - view exists: v2_attachment_site_candidates: observed=1; expected=1
* **PASS** - view exists: v2_attachment_site_high_priority: observed=1; expected=1
* **PASS** - view exists: v2_attachment_site_summary: observed=1; expected=1
* **PASS** - view exists: Virus_Proteins: observed=1; expected=1
* **PASS** - view exists: Ligand_Atoms_Smiles: observed=1; expected=1
* **PASS** - view exists: Functional_GROUPED: observed=1; expected=1
* **PASS** - view exists: ligand_atoms: observed=1; expected=1
* **PASS** - view exists: solvent_exposed_atoms: observed=1; expected=1
* **PASS** - view exists: RUPLEY_SASA_DATA: observed=1; expected=1
* **PASS** - view exists: SMILES_MAP_PDB: observed=1; expected=1
* **PASS** - view exists: Functional_Group_Atoms: observed=1; expected=1
* **PASS** - view exists: Arpeggio_Contacts_Data: observed=1; expected=1
* **PASS** - view exists: receptor_binding_pocket: observed=1; expected=1
* **PASS** - view exists: Covalent_Noncovalent: observed=1; expected=1
* **PASS** - view exists: distal_atoms: observed=1; expected=1
* **PASS** - all expected compatibility views present: observed=23; expected=23

## Code freeze

* **PASS** - 07_map_cif_atoms_to_smiles.py VERSION: observed=legacy_mcs_etkdg_uff_cif_v2.5; expected=legacy_mcs_etkdg_uff_cif_v2.5
* **PASS** - 08_calculate_ligand_sasa.py VERSION: observed=biopython-shrake_rupley-1.40-cif-v2.1; expected=biopython-shrake_rupley-1.40-cif-v2.1
* **PASS** - 09_run_arpeggio.py VERSION: observed=arpeggio-cif-v2.2; expected=arpeggio-cif-v2.2
* **PASS** - 10_calculate_ligand_geometry.py VERSION: observed=cif-ligand-geometry-v2.4; expected=cif-ligand-geometry-v2.4
* **PASS** - 11_assign_functional_groups.py VERSION: observed=rdkit-smarts-functional-groups-v2.3; expected=rdkit-smarts-functional-groups-v2.3
* **PASS** - 12_build_protacability.py VERSION: observed=protacability-cif-v2.8; expected=protacability-cif-v2.8
* **PASS** - 13_build_attachment_sites.py VERSION: observed=attachment-sites-cif-v2.6; expected=attachment-sites-cif-v2.6
* **PASS** - 14_build_compatibility_views.py VERSION: observed=compatibility-views-cif-v2.7; expected=compatibility-views-cif-v2.7

## Foundation

* **PASS** - frozen mmCIF manifest rows: observed=11533; expected=11533
* **PASS** - structure rows: observed=7610; expected=7610
* **PASS** - unique PDB entries: observed=7610; expected=7610
* **PASS** - retained ligand instances: observed=7355; expected=7355
* **PASS** - resolved-chemistry denominator: observed=7335; expected=7335
* **PASS** - included instances without resolved chemistry: observed=20; expected=20

## Stage 07 mapping

* **PASS** - mapping population: observed=7335; expected=7335
* **PASS** - mapping status distribution: observed={'complete': 3864, 'complete_altloc_resolved': 570, 'partial_ccd_difference': 2182, 'skipped_pending_remediation': 719}; expected={'complete': 3864, 'complete_altloc_resolved': 570, 'partial_ccd_difference': 2182, 'skipped_pending_remediation': 719}
* **PASS** - downstream-eligible mappings: observed=6616; expected=6616
* **PASS** - pending remediation population: observed=719; expected=719
* **PASS** - failed/time-out mappings: observed=0; expected=0
* **PASS** - unique SMILES index per occurrence: observed=0; expected=0

## Stage 08 SASA

* **PASS** - instance coverage: observed=7355; expected=7355
* **PASS** - 1.40 A / 100-point / HOH-removed parameter consistency: observed=0; expected=0
* **PASS** - selected ligand atoms missing SASA: observed=0; expected=0

## Stage 09 Arpeggio

* **PASS** - latest completed outcomes in included release population: observed=7355; expected=7355
* **PASS** - latest unresolved outcomes in included release population: observed=0; expected=0
* **PASS** - duplicate latest outcome keys in included release population: observed=0; expected=0
* **PASS** - included instances without Arpeggio outcome: observed=0; expected=0

## Stage 10 geometry

* **PASS** - inventory coverage: observed=7355; expected=7355
* **PASS** - protein-applicable instances: observed=7308; expected=7308
* **PASS** - no-protein not-applicable instances: observed=47; expected=47
* **PASS** - unexpected geometry statuses: observed=0; expected=0

## Stage 11 functional groups

* **PASS** - complete summaries: observed=7335; expected=7335
* **PASS** - SMARTS match count: observed=53224; expected=53224
* **PASS** - mapped functional-group atom occurrences: observed=134419; expected=134419
* **PASS** - unmapped functional-group atom occurrences: observed=24472; expected=24472
* **PASS** - unmapped FG atoms on complete/altloc-complete mappings: observed=0; expected=0
* **PASS** - mapped FG atoms lacking element validation: observed=0; expected=0

## Stage 12 PROTACability

* **PASS** - warhead/linkability coverage: observed=7355; expected=7355
* **PASS** - target-context distribution: observed={'applicable_contacting_protein_chain': 6784, 'not_applicable_no_contacting_protein_chain': 524, 'not_applicable_no_protein_atoms': 47}; expected={'applicable_contacting_protein_chain': 6784, 'not_applicable_no_contacting_protein_chain': 524, 'not_applicable_no_protein_atoms': 47}
* **PASS** - assessed ligand instances: observed=6784; expected=6784
* **PASS** - target-chain assessment rows: observed=9462; expected=9462
* **PASS** - readiness instance coverage: observed=6784; expected=6784
* **PASS** - readiness row count: observed=9462; expected=9462
* **PASS** - readiness tier distribution: observed={'Exploratory degrader-design readiness': 173, 'High degrader-design readiness': 6579, 'Moderate degrader-design readiness': 2289, 'Weak degrader-design readiness': 421}; expected={'High degrader-design readiness': 6579, 'Moderate degrader-design readiness': 2289, 'Weak degrader-design readiness': 421, 'Exploratory degrader-design readiness': 173}
* **PASS** - all target rows supported by direct Stage-09 protein contact: observed=0; expected=0
* **PASS** - surface-lysine score recomputes from NZ coverage/SASA only: observed=0; expected=0
* **PASS** - readiness score recomputes from warhead + surface lysines only: observed=0; expected=0
* **PASS** - ligand-to-lysine proximity disabled in readiness rows: observed=0; expected=0
* **PASS** - zero accessible NZ cannot be Moderate/High: observed=0; expected=0
* **PASS** - exactly one best target chain per assessed ligand: observed=0; expected=0

## Stage 13 attachment sites

* **PASS** - complete instance summaries: observed=7355; expected=7355
* **PASS** - atom rows: observed=227080; expected=227080
* **PASS** - candidate atoms: observed=62392; expected=62392
* **PASS** - High-priority atoms: observed=1544; expected=1544
* **PASS** - instances with candidates: observed=5814; expected=5814
* **PASS** - instances with High sites: observed=1120; expected=1120
* **PASS** - priority tier distribution: observed={'Exploratory attachment-site priority': 53606, 'High attachment-site priority': 1544, 'Low attachment-site priority': 164688, 'Moderate attachment-site priority': 7242}; expected={'Low attachment-site priority': 164688, 'Exploratory attachment-site priority': 53606, 'Moderate attachment-site priority': 7242, 'High attachment-site priority': 1544}
* **PASS** - atom chemical-role distribution: observed={'conditional_substitution_site': 18530, 'direct_attachment_atom': 9072, 'functional_group_context_only': 83144, 'unclassified_atom_context': 116334}; expected={'unclassified_atom_context': 116334, 'functional_group_context_only': 83144, 'conditional_substitution_site': 18530, 'direct_attachment_atom': 9072}
* **PASS** - candidate-core rule equivalence: observed=0; expected=0
* **PASS** - High-rule violations: observed=0; expected=0
* **PASS** - High tier/flag consistency: observed=0; expected=0

## Occurrence traceability

* **PASS** - cross-instance mapping atom joins: observed=0; expected=0
* **PASS** - cross-instance SASA atom joins: observed=0; expected=0
* **PASS** - cross-instance geometry atom joins: observed=0; expected=0
* **PASS** - cross-instance functional-group atom joins: observed=0; expected=0
* **PASS** - cross-instance attachment-site atom joins: observed=0; expected=0
* **PASS** - cross-instance Arpeggio ligand-atom joins: observed=0; expected=0

## Stage 14 compatibility

* **PASS** - latest build completed 23/23/0: observed=('completed', 23, 23, 0); expected=('completed', 23, 23, 0)
* **PASS** - v2_target_lysine_accessibility pinned to current evidence: observed=all tokens present; expected=all tokens present
* **PASS** - v2_protacability_target_context pinned to current evidence: observed=all tokens present; expected=all tokens present
* **PASS** - v2_protacability_best pinned to current evidence: observed=all tokens present; expected=all tokens present
* **PASS** - v2_attachment_site_candidates pinned to current evidence: observed=all tokens present; expected=all tokens present
* **PASS** - v2_attachment_site_high_priority pinned to current evidence: observed=all tokens present; expected=all tokens present
* **PASS** - v2_attachment_site_summary pinned to current evidence: observed=all tokens present; expected=all tokens present

## 3EKY/DR7 regression

* **PASS** - unique included DR7 occurrence: observed=1; expected=1
* **PASS** - 51/51 authoritative mapping: observed=('complete', 51, 51, 51, 1.0, 1); expected=('complete', 51, 51, 51, 1.0, 1)
* **PASS** - functional-group mapped atoms element-validated: observed=0; expected=0
* **PASS** - direct target chains A;B: observed={'target_context_status': 'applicable_contacting_protein_chain', 'contacting_protein_chain_count': 2, 'contacting_protein_chain_ids': 'A;B', 'target_chain_selection_basis': 'stage09_arpeggio_direct_protein_contact'}; expected={'status': 'applicable_contacting_protein_chain', 'count': 2, 'chains': 'A;B'}
* **PASS** - surface-lysine scores A=52.86, B=81.43: observed={'A': 52.86, 'B': 81.43}; expected={'A': 52.86, 'B': 81.43}
* **PASS** - no forced High direct-attachment site: observed=0; expected=0
* **PASS** - conditional Moderate sites: observed=['CAO', 'CAR', 'CAS', 'NBD']; expected=['CAO', 'CAR', 'CAS', 'NBD']

## Database health

* **PASS** - PRAGMA integrity_check: observed=ok; expected=ok
* **PASS** - PRAGMA foreign_key_check: observed=0; expected=0

## Frozen release snapshot for manuscript / reviewer response

* Frozen mmCIF files: 11533
* Unique PDB entries: 7610
* Retained ligand instances: 7355
* Resolved-chemistry instances: 7335
* Stage-07 downstream-usable mappings: 6616
* Stage-07 pending remediation instances: 719
* Stage-09 completed Arpeggio outcomes in included release population: 7355
* Excluded/non-release ligand instances with completed Arpeggio history: 4 (reported for provenance; outside the 7,355 release denominator)
* Stage-10 protein-applicable geometry: 7308
* Stage-10 no-protein target N/A: 47
* Stage-12 directly contacting target contexts: 6784
* Stage-12 target-chain assessments: 9462
* Stage-13 ligand instances with candidate attachment sites: 5814
* Stage-13 ligand instances with High direct attachment sites: 1120
* Stage-13 candidate atoms: 62392
* Stage-13 High-priority direct attachment atoms: 1544

### Interpretation guardrails

* The 7,355 retained-ligand denominator covers the included structure-derived release population; Arpeggio history for excluded/non-release occurrences is retained for provenance but is not counted in that denominator. The 7,335 denominator is the resolved-chemistry population eligible for SMILES-based mapping and SMARTS annotation.
* Target PROTACability uses Stage-09-confirmed ligand-contacting protein chains, then evaluates lysine NZ solvent accessibility across the entire selected chain surface.
* Ligand-to-lysine distance is retained only in upstream descriptive geometry where present and is not used in the v2.8 target-accessibility or degrader-readiness score.
* High attachment-site priority requires a mapped, solvent-exposed, low-strong-contact atom with direct atom-level attachment chemistry, an outward-facing vector, and a locally clear corridor.
* These are structural-priority heuristics for follow-up design, not experimentally calibrated predictors of linker tolerance, ubiquitination, or degradation.

