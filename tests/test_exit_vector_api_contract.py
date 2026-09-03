import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "RANDY"))

import vlismod_data_routes as routes  # noqa: E402


class ExitVectorApiContractTests(unittest.TestCase):
    def test_exit_vector_fields_survive_row_merge_and_target_snapshot(self):
        assessment = {
            "virus_name": "HIV_1",
            "protein_type": "protease",
            "pdb_code": "3EKY",
            "chain_id": "A",
            "candidate_ligand_resnames": "DR7",
        }
        readiness = {
            **assessment,
            "ligand_exit_geometry_score": 53.33,
            "clear_exit_candidate_count": 6,
            "best_linker_geometry_class": "Open outward ligand exit-vector cue",
        }
        warhead = {
            "pdb_code": "3EKY",
            "ligand_resname": "DR7",
            "ligand_chain": "A",
            "ligand_residue_id": "100",
            "candidate_linker_atom_count": 12,
            "outward_supported_candidate_count": 7,
            "clear_exit_candidate_count": 6,
            "ligand_exit_geometry_score": 53.33,
        }

        merged = routes._merge_optional_protacability_data(
            [assessment], readiness_rows=[readiness], warhead_rows=[warhead]
        )[0]
        self.assertEqual(merged["ligand_exit_geometry_score"], 53.33)
        self.assertEqual(merged["clear_exit_candidate_count"], 6)
        self.assertEqual(merged["outward_supported_candidate_count"], 7)
        self.assertEqual(merged["candidate_linker_atom_count"], 12)
        self.assertEqual(
            merged["best_linker_geometry_class"], "Open outward ligand exit-vector cue"
        )

        snapshot = routes._protacability_enrichment_snapshot(merged)
        self.assertEqual(snapshot["ligand_exit_geometry_score"], 53.33)
        self.assertEqual(snapshot["clear_exit_candidate_count"], 6)
        self.assertEqual(snapshot["outward_supported_candidate_count"], 7)


if __name__ == "__main__":
    unittest.main()
