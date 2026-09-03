import unittest
from pathlib import Path


class TargetDetailStructureViewerTests(unittest.TestCase):
    def test_target_structure_rows_load_their_own_viewer_payload(self):
        template = Path("templates/protacability_assessment.html").read_text()
        self.assertIn('class="protac-target-structure-row', template)
        self.assertIn('data-pdb-code=', template)
        self.assertIn("const loadStructureModel = function(row)", template)
        self.assertIn('attachProtacViewerControls(pdbCode, chainId, ligandResname, ligandChain, ligandResidueId, ligandResname, candidateLinkerAtomIds, {}, syncSelectedLigand)', template)
        self.assertIn('id="protac-selected-context"', template)
        self.assertIn('api/protacability/structure_detail/', template)
        self.assertIn('renderLigandContextSelector(detail.ligand_contexts || [], selectedLigand)', template)

    def test_target_structure_rows_are_keyboard_accessible(self):
        template = Path("templates/protacability_assessment.html").read_text()
        self.assertIn('role="button" tabindex="0"', template)
        self.assertIn("event.key === 'Enter' || event.key === ' '", template)

    def test_start_analysis_tracks_the_selected_ligand_occurrence(self):
        template = Path("templates/protacability_assessment.html").read_text()
        self.assertIn('id="protac-start-analysis"', template)
        self.assertIn("const updateAnalysisContext = function(selection = {})", template)
        self.assertIn("'data-analysis-pdb-code': String(selection.pdbCode || activePdbCode || '')", template)
        self.assertIn("'data-analysis-ligand-instance-id': String(selection.ligandInstanceId || activeLigandInstanceId || '')", template)
        self.assertIn('updateAnalysisContext(selection);', template)


if __name__ == '__main__':
    unittest.main()
