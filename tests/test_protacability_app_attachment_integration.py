import sqlite3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "RANDY"))
sys.path.insert(0, str(ROOT / "TOOLS"))

import vlismod_data_routes as routes  # noqa: E402
from protacability_attachment_schema import CREATE_STATEMENTS  # noqa: E402


class QueryArgs(dict):
    def get(self, key, default=None, type=None):
        value = super().get(key, default)
        if type is not None and value not in (None, ""):
            return type(value)
        return value

    def getlist(self, key):
        value = self.get(key)
        if value is None:
            return []
        return value if isinstance(value, list) else [value]


class ProtacabilityAttachmentAppIntegrationTests(unittest.TestCase):
    def test_decorated_rows_include_exact_attachment_summary_and_filter(self):
        assessment_rows = [
            {
                "virus_name": "Test virus",
                "protein_type": "Protease",
                "pdb_code": "3EKY",
                "chain_id": "A",
                "candidate_ligand_resnames": "DR7",
                "has_candidate_ligand": 1,
                "has_exposed_lysine": 1,
                "has_ligand_proximal_exposed_lysine": 1,
                "protacability_proxy_score": 80,
                "protacability_tier": "High structural priority",
                "exposed_lys_count": 3,
            }
        ]
        warhead_rows = [
            {
                "virus_name": "Test virus",
                "protein_type": "Protease",
                "pdb_code": "3EKY",
                "model_id": 0,
                "ligand_resname": "DR7",
                "ligand_chain": "A",
                "ligand_residue_id": 101,
                "ligand_insertion_code": "",
                "candidate_linker_atom_count": 2,
                "candidate_linker_atom_ids": "1567;1568",
                "solvent_exposed_ligand_atom_count": 4,
                "pdb_to_smiles_mapped_atom_count": 12,
                "rdkit_valid_smiles": 1,
            }
        ]
        attachment_rows = [
            {
                "analysis_id": 42,
                "pdb_code": "3EKY",
                "model_id": 0,
                "ligand_chain": "A",
                "ligand_residue_id": 101,
                "ligand_insertion_code": "",
                "ligand_resname": "DR7",
                "analysis_status": "completed",
                "mapping_status": "complete",
                "eligibility_status": "full_analysis_eligible",
                "has_attachment_site_evidence": 1,
                "attachment_region_count": 2,
                "candidate_atom_count": 5,
                "best_attachment_score": 0.87,
                "best_attachment_confidence": "high",
                "method_version": routes.ATTACHMENT_METHOD_VERSION,
                "instance_resolution_status": "resolved",
                "instance_ambiguity_flag": 0,
            }
        ]

        rows = routes._decorate_protacability_rows(
            assessment_rows,
            readiness_rows=[],
            warhead_rows=warhead_rows,
            attachment_rows=attachment_rows,
        )

        self.assertEqual(rows[0]["attachment_analysis_id"], 42)
        self.assertEqual(rows[0]["attachment_region_count"], 2)
        self.assertEqual(rows[0]["attachment_candidate_atom_count"], 5)
        self.assertEqual(rows[0]["has_attachment_site_evidence"], 1)

        filters = routes._build_protacability_filters(QueryArgs({"has_attachment_sites": "1"}))
        self.assertEqual(len(routes._filter_protacability_rows(rows, filters)), 1)
        filters = routes._build_protacability_filters(QueryArgs({"has_attachment_sites": "0"}))
        self.assertEqual(len(routes._filter_protacability_rows(rows, filters)), 0)

    def test_attachment_detail_payload_serializes_regions_and_candidate_serials(self):
        conn = sqlite3.connect(":memory:")
        conn.row_factory = sqlite3.Row
        for statement in CREATE_STATEMENTS:
            conn.execute(statement)
        conn.execute(
            """
            INSERT INTO protacability_attachment_analysis (
                analysis_id, pdb_code, model_id, ligand_chain, ligand_residue_id,
                ligand_insertion_code, ligand_resname, graph_id, analysis_status,
                mapping_status, source_data_completeness_json, has_attachment_site_evidence,
                attachment_region_count, candidate_atom_count, best_attachment_score,
                best_attachment_confidence, bond_source, sasa_source, contact_source,
                eligibility_status, software_versions_json, method_version,
                calculation_parameters_json, generated_at
            )
            VALUES (
                7, '3EKY', 0, 'A', 101, '', 'DR7', NULL, 'completed',
                'complete', '{}', 1, 1, 2, 0.91, 'high', 'graph',
                'sasa', 'contacts', 'full_analysis_eligible', '{}', ?, '{}',
                '2026-07-31T00:00:00Z'
            )
            """,
            (routes.ATTACHMENT_METHOD_VERSION,),
        )
        conn.execute(
            """
            INSERT INTO protacability_attachment_regions (
                analysis_id, region_id, member_atom_ids_json, member_smiles_indices_json,
                candidate_atom_ids_json, best_candidate_atom_id, interaction_summary_json,
                region_score, confidence, reasons_json, cautions_json, method_version
            )
            VALUES (7, 'R1', '[1567,1568]', '[3,4]', '[1567]', 1567, '{}', 0.9, 'high', '["exposed"]', '[]', ?)
            """,
            (routes.ATTACHMENT_METHOD_VERSION,),
        )
        conn.execute(
            """
            INSERT INTO protacability_attachment_atoms (
                analysis_id, region_id, pdb_atom_serial, pdb_atom_name, element,
                smiles_atom_index, interaction_types_json, functional_group_annotations_json,
                candidate_attachment_flag, surface_defining_flag, attachment_score,
                confidence, reasons_json, cautions_json, method_version
            )
            VALUES
                (7, 'R1', 1567, 'CAO', 'C', 3, '[]', '[]', 1, 1, 0.91, 'high', '["terminal"]', '[]', ?),
                (7, 'R1', 1568, 'CAS', 'S', 4, '[]', '[]', 1, 0, 0.89, 'high', '["exposed"]', '[]', ?)
            """,
            (routes.ATTACHMENT_METHOD_VERSION, routes.ATTACHMENT_METHOD_VERSION),
        )

        payload = routes._attachment_detail_payload(
            conn,
            {
                "pdb_code": "3EKY",
                "model_id": 0,
                "ligand_chain": "A",
                "ligand_residue_id": 101,
                "ligand_insertion_code": "",
                "ligand_resname": "DR7",
            },
        )

        self.assertTrue(payload["data_available"])
        self.assertEqual(payload["summary"]["attachment_analysis_id"], 7)
        self.assertEqual(payload["candidate_atom_serials"], [1567, 1568])
        self.assertEqual(payload["surface_atom_serials"], [1567])
        self.assertEqual(payload["regions"][0]["reasons"], ["exposed"])
        self.assertEqual(payload["regions"][0]["candidate_atom_serials"], [1567, 1568])
        self.assertEqual(payload["regions"][0]["surface_atom_serials"], [1567])
        self.assertEqual(payload["atoms"][0]["pdb_atom_name"], "CAO")


if __name__ == "__main__":
    unittest.main()
