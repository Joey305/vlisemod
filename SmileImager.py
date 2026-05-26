from rdkit import Chem
from rdkit.Chem import Draw

# Define the SMILES string
smiles = "COC(=O)N[C@H](C(=O)NN(Cc1ccc(-c2ccccn2)cc1)C[C@H](O)[C@H](Cc1ccccc1)NC(=O)[C@@H](NC(=O)COCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCNC(=O)CC(C)(C)CC(=O)N1C[C@H](O)C[C@H]1C(=O)NCc1ccc(-c2cnco2)cc1)C(C)(C)C)C(C)(C)C"

# Convert SMILES to a molecule object
mol = Chem.MolFromSmiles(smiles)

# Draw molecule with a transparent background in SVG format
drawer = Draw.MolDraw2DSVG(1000, 1000)  # High resolution
drawer.DrawMolecule(mol)
drawer.FinishDrawing()

# Save as SVG file
with open("molecule.svg", "w") as f:
    f.write(drawer.GetDrawingText())

print("SVG file saved as 'molecule.svg'")
