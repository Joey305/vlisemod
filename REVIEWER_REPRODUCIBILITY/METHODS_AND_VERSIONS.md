# Canonical methods and release generations

| Stage | Script | Frozen release generation |
| --- | --- | --- |
| 07 | `07_map_cif_atoms_to_smiles.py` | `legacy_mcs_etkdg_uff_cif_v2.5` |
| 08 | `08_calculate_ligand_sasa.py` | `biopython-shrake_rupley-1.40-cif-v2.1` |
| 09 | `09_run_arpeggio.py` | `arpeggio-cif-v2.2` |
| 10 | `10_calculate_ligand_geometry.py` | `cif-ligand-geometry-v2.4` |
| 11 | `11_assign_functional_groups.py` | `rdkit-smarts-functional-groups-v2.3` |
| 12 | `12_build_protacability.py` | `protacability-cif-v2.8` |
| 13 | `13_build_attachment_sites.py` | `attachment-sites-cif-v2.6` |
| 14 | `14_build_compatibility_views.py` | `compatibility-views-cif-v2.7` |
| 15 | `15_validate_database.py` | `final-release-validation-cif-v2.4` |

Stages 01–03 inventory, create, and ingest the frozen CIF corpus. Stage 05 applies retained-ligand curation. Stage 06 loads the frozen component-chemistry table; Stage 06a loads the explicit pending-remediation registry. The later stages are then executed exactly in numerical order.

The historical `04_identify_ligand_instances.py` is not a separate scientific stage in this package because it is a compatibility wrapper for the Stage-03 ingest logic. The following noncanonical scripts are deliberately absent: `10_calculate_ligand_geometry_OLD.py` and `13_build_attachment_site.py`.

Stage 15 validates the complete release target: 11,533 CIF files; 7,610 structures; 7,355 retained ligand instances; 7,335 resolved chemistry instances; mapping counts 3,864 complete / 570 altloc complete / 2,182 partial / 719 pending remediation; 7,355 completed Arpeggio runs; 7,308 geometry-complete cases; 53,224 functional-group matches; 9,462 target-chain assessment rows; 227,080 attachment-atom rows; and 23 compatibility views. It is intentionally not weakened for partial runs.
