import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "TOOLS" / "assess_exit_vector_cue.py"
SPEC = importlib.util.spec_from_file_location("assess_exit_vector_cue", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def atom(**overrides):
    base = {
        "candidate_attachment_atom": 1,
        "points_away_from_pocket": 1,
        "local_corridor_clear": 1,
        "forward_clearance_a": 8.0,
        "attachment_priority_score": 85.0,
        "atom_site_id": "1",
        "exact_atom": "C1",
        "element": "C",
        "smiles_atom_indices": "0",
        "attachment_priority_tier": "High attachment-site priority",
        "chemical_context": "direct_handle_context",
        "solvent_exposed": 1,
        "strong_contact_count": 0,
        "outward_score": 0.5,
        "forward_obstruction_count": 0,
    }
    base.update(overrides)
    return base


class ExitVectorCueTests(unittest.TestCase):
    def test_clear_outward_candidates_produce_the_recommended_cue(self):
        result = MODULE.assess_atoms([atom(), atom(atom_site_id="2", local_corridor_clear=0)])
        self.assertEqual(result["candidate_attachment_atom_count"], 2)
        self.assertEqual(result["outward_candidate_count"], 2)
        self.assertEqual(result["locally_clear_outward_candidate_count"], 1)
        self.assertEqual(result["exit_vector_score"], 70.0)
        self.assertEqual(result["cue"], "Locally clear ligand exit-vector cue")

    def test_outward_but_obstructed_candidates_are_not_called_clear(self):
        result = MODULE.assess_atoms([atom(local_corridor_clear=0)])
        self.assertEqual(result["exit_vector_score"], 40.0)
        self.assertEqual(result["cue"], "Outward ligand atom cue with local forward obstruction")

    def test_no_candidates_has_zero_score(self):
        result = MODULE.assess_atoms([atom(candidate_attachment_atom=0)])
        self.assertEqual(result["exit_vector_score"], 0.0)
        self.assertEqual(result["cue"], "No candidate ligand attachment atom")
