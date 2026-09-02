import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("vlismod_pymol_identity_app", ROOT / "app.py")
vlismod_app = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = vlismod_app
SPEC.loader.exec_module(vlismod_app)


class PymolOccurrenceIdentityTests(unittest.TestCase):
    def test_rebuilt_sasa_and_protein_atoms_use_stable_names_not_cif_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            previous = Path.cwd()
            try:
                # The writer intentionally uses a relative output directory.
                import os
                os.chdir(tmp)
                path = vlismod_app.write_pymol_script(
                    "3EKY", "DR7", "A", "100",
                    {
                        "solvent_exposed_atoms": [{"atom_id": "1517", "chain": "A", "exact_atom": "CAA"}],
                        "binding_pocket": [{"residue_chain": "B", "residue_number": "42", "residue_atom": "CA"}],
                    },
                )
                script = Path(path).read_text()
            finally:
                os.chdir(previous)

        self.assertIn("resn DR7 and chain A and resi 100 and name CAA", script)
        self.assertIn("polymer and chain B and resi 42 and name CA", script)
        self.assertNotIn("chain A and id 1517", script)

    def test_index_form_posts_the_selected_occurrence_id(self):
        form = (ROOT / "templates" / "index.html").read_text()
        scripts = (ROOT / "static" / "js" / "scripts.js").read_text()
        self.assertIn('name="ligand_instance_id"', form)
        self.assertIn("dataset.ligandInstanceId", scripts)
        self.assertIn("$('#ligand_instance_id').val", scripts)


if __name__ == "__main__":
    unittest.main()
