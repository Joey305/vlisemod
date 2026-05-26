from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import os
import csv
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO)

# Function to parse PDB file and extract atom data while preserving atom IDs
def parse_pdb_file(pdb_file, ligand_code):
    # Convert PosixPath to string
    pdb_file_str = str(pdb_file)
    
    mol = Chem.MolFromPDBFile(pdb_file_str, removeHs=False)
    if mol is None:
        logging.error(f"Failed to parse PDB file {pdb_file} for ligand {ligand_code}.")
        return None, None

    # Extract atoms specific to the ligand while preserving PDB atom IDs
    ligand_atoms = []
    atom_mapping = {}  # Maps RDKit atom index to PDB atom ID
    for atom in mol.GetAtoms():
        res_name = atom.GetPDBResidueInfo().GetResidueName().strip()
        if res_name == ligand_code:
            pdb_atom_id = atom.GetPDBResidueInfo().GetSerialNumber()
            atom_idx = atom.GetIdx()
            ligand_atoms.append(atom)
            atom_mapping[atom_idx] = pdb_atom_id

    logging.info(f"Extracted {len(ligand_atoms)} atoms from PDB file {pdb_file} for ligand {ligand_code}.")
    return mol, atom_mapping

# Function to generate a 3D conformer from a SMILES string
def generate_conformer(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens to the molecule
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # Generate a 3D conformer
    AllChem.UFFOptimizeMolecule(mol)  # Optimize the conformer using UFF
    logging.info("Conformer generated successfully.")
    return mol

# Function to match atoms between PDB and SMILES using bond connectivity
def match_atoms_by_bond_connectivity(pdb_mol, smiles_mol, atom_mapping):
    # Use RDKit's MCS (Maximum Common Substructure) to find the best match
    res = rdFMCS.FindMCS([pdb_mol, smiles_mol], completeRingsOnly=True, bondCompare=rdFMCS.BondCompare.CompareOrder)
    patt = Chem.MolFromSmarts(res.smartsString)
    pdb_matches = pdb_mol.GetSubstructMatch(patt)
    smiles_matches = smiles_mol.GetSubstructMatch(patt)

    if len(pdb_matches) == 0 or len(smiles_matches) == 0:
        logging.info("No matches found between PDB and SMILES conformer.")
        return []

    # Use the atom_mapping to retain the correct PDB atom IDs
    matches = [(atom_mapping[pdb_idx], pdb_idx, smiles_idx) for pdb_idx, smiles_idx in zip(pdb_matches, smiles_matches)]
    logging.info(f"Found {len(matches)} matches between PDB atoms and SMILES conformer atoms based on bond connectivity.")
    return matches

# Function to save the matched atoms to a CSV file
def save_matched_atoms_to_csv(matches, pdb_id, output_csv_file):
    with open(output_csv_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        for match in matches:
            writer.writerow([pdb_id] + list(match))
    logging.info(f"Matches saved to {output_csv_file}")

# Function to process the input .txt file and map atoms
def process_txt_input(base_input_dir, output_csv_file, virus_names):
    for virus_name in virus_names:
        # Construct the paths for the TXT file and PDB directory
        txt_file = Path(base_input_dir) / 'Sorted_PDB_Files' / virus_name / 'with_ligands_smiles_functional_groups.txt'
        pdb_dir = Path(base_input_dir) / virus_name 

        if not txt_file.exists():
            logging.error(f"TXT file not found: {txt_file}")
            continue

        with open(txt_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)

            for row in reader:
                if len(row) < 3:
                    logging.error(f"Incorrect format in row: {row}")
                    continue

                pdb_id = row[0]
                ligand_code = row[1]
                smiles = row[2]

                logging.info(f"Processing PDB ID: {pdb_id} with Ligand: {ligand_code}")

                pdb_file = pdb_dir / f"{pdb_id}.pdb"
                if not pdb_file.exists():
                    logging.error(f"PDB file not found: {pdb_file}")
                    continue

                pdb_mol, atom_mapping = parse_pdb_file(pdb_file, ligand_code)
                if pdb_mol is None:
                    continue

                conformer = generate_conformer(smiles)
                matches = match_atoms_by_bond_connectivity(pdb_mol, conformer, atom_mapping)

                # Save matches to CSV
                save_matched_atoms_to_csv(matches, pdb_id, output_csv_file)

# Main function to orchestrate the processing
def main():
    base_input_dir = './Database_DATA'  # Base directory containing virus subdirectories
    output_csv_file = 'matched_atoms_output.csv'

    virus_names = [
        'Human immunodeficiency virus 1',
        # Add more virus names as needed
    ]

    # Create the CSV file with headers
    with open(output_csv_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['PDB ID', 'PDB Atom ID', 'PDB Atom Index', 'SMILES Atom Index'])

    # Process the input TXT file
    process_txt_input(base_input_dir, output_csv_file, virus_names)

if __name__ == "__main__":
    main()
