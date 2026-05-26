import os
import csv
from rdkit import Chem
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    filename='functional_groups_log.log',
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.DEBUG
)

# Function to load functional groups from a text file
def load_functional_groups(file_path):
    functional_groups_smarts = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip():  # Skip empty lines
                group, smarts = line.strip().split(':')
                functional_groups_smarts[group.strip()] = smarts.strip()
    return functional_groups_smarts

# Function to identify functional groups in a molecule and match them to atom indices
def identify_functional_groups(smiles, functional_groups_smarts):
    if smiles.lower() == 'n/a':
        logging.error(f"Skipping invalid SMILES: {smiles}")
        return []

    try:
        mol = Chem.MolFromSmiles(smiles)
    except ValueError as e:
        logging.error(f"Error processing SMILES: {smiles} - {str(e)}")
        return []

    if mol is None:
        logging.error(f"Failed to parse SMILES: {smiles}")
        return []

    matched_groups = []
    for group_name, smarts in functional_groups_smarts.items():
        fg_mol = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(fg_mol):
            matches = mol.GetSubstructMatches(fg_mol)
            for match in matches:
                matched_groups.append((group_name, match))
    
    return matched_groups

# Helper function to convert integer to Roman numeral
def int_to_roman(n):
    roman_numerals = {
        1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V',
        6: 'VI', 7: 'VII', 8: 'VIII', 9: 'IX', 10: 'X'
    }
    return roman_numerals.get(n, str(n))

# Function to process the SMILES .txt file and create a new text file with matched functional groups and their atom IDs
def process_smiles_file(input_file, output_file, matched_atoms_file, functional_groups_smarts):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(matched_atoms_file, 'r') as matched_file:
        reader = csv.reader(infile, delimiter='\t')
        matched_reader = csv.reader(matched_file)

        # Skip the header row in both the input and matched atoms files
        next(reader)
        next(matched_reader)
        
        # Load the matched atom index to PDB atom ID, chain, exact atom, and atom type mappings
        matched_atoms = {
            (row[0], int(row[4])): (row[1], row[2], row[5], row[6]) for row in matched_reader
        }

        # Write the header to the output text file
        header = "PDB ID\tLigand ID\tChain\tFunctional Group\tAtom ID\tExact Atom\tAtom Type\n"
        outfile.write(header)

        # Dictionary to keep track of how many times a functional group has been encountered
        group_counter = {}

        for row in reader:
            pdb_id = row[0]  # Assuming PDB ID is in the first column
            ligand_id = row[1]  # Assuming Ligand ID is in the second column
            smiles = row[-1]  # Assuming SMILES is in the last column

            matched_groups = identify_functional_groups(smiles, functional_groups_smarts)
            
            if not matched_groups:
                logging.warning(f"No functional groups matched for SMILES: {smiles} (PDB ID: {pdb_id}, Ligand ID: {ligand_id})")
                continue

            for group_name, atom_indices in matched_groups:
                group_key = (pdb_id, ligand_id, group_name)
                count = group_counter.get(group_key, 0) + 1
                group_counter[group_key] = count
                distinguished_group_name = f"{group_name} ({int_to_roman(count)})"

                for idx in atom_indices:
                    if (pdb_id, idx) in matched_atoms:
                        atom_id, chain, exact_atom, atom_type = matched_atoms[(pdb_id, idx)]
                        # Write each row to the text file
                        outfile.write(f"{pdb_id}\t{ligand_id}\t{chain}\t{distinguished_group_name}\t{atom_id}\t{exact_atom}\t{atom_type}\n")
                        logging.info(
                            f"Processed functional group '{distinguished_group_name}' for ligand with SMILES: {smiles}, "
                            f"PDB ID: {pdb_id}, Ligand ID: {ligand_id}, Chain: {chain}, Atom ID: {atom_id}, Exact Atom: {exact_atom}, Atom Type: {atom_type}"
                        )

    logging.info(f"Processing complete. Output saved to {output_file}")

# Function to process all SMILES files based on virus names
def process_virus_smiles_files(base_input_dir, virus_names, functional_groups_smarts, matched_atoms_file):
    for virus_name in virus_names:
        smiles_file = Path(base_input_dir) / 'Sorted_PDB_Files' / virus_name / 'with_ligands_smiles.txt'
        output_file = smiles_file.with_name(smiles_file.stem + '_functional_groups.txt')
        
        if not smiles_file.exists():
            logging.error(f"SMILES file not found for virus '{virus_name}': {smiles_file}")
            continue

        logging.info(f"Processing SMILES file for virus '{virus_name}': {smiles_file}")
        process_smiles_file(str(smiles_file), str(output_file), matched_atoms_file, functional_groups_smarts)

# Main function
def main():
    base_input_dir = './Database_DATA'
    functional_groups_file = 'functional_groups.txt'  # Path to your functional groups file
    matched_atoms_file = 'matched_atoms_output.csv'  # Path to your matched atoms file

    virus_names = [
        'Human immunodeficiency virus 1',
        "Severe acute respiratory syndrome coronavirus 2",
        "Human papillomavirus",
        "Human papillomavirus 1",
        "Human papillomavirus 11",
        "Human papillomavirus 16",
        "Human papillomavirus 18",
        "Human papillomavirus 26",
        "Human papillomavirus 31",
        "Human papillomavirus 33",
        "Human papillomavirus 35",
        "Human papillomavirus 4",
        "Human papillomavirus 45",
        "Human papillomavirus 49",
        "Human papillomavirus 51",
        "Human papillomavirus 52",
        "Human papillomavirus 53",
        "Human papillomavirus 58",
        "Human papillomavirus 59",
        "Human papillomavirus 6",
        "Human papillomavirus 66",
        "Human papillomavirus type 16",
        "Human papillomavirus type 6",
        "Human papillomavirus type 6a"
    ]

    # Load functional groups from file
    functional_groups_smarts = load_functional_groups(functional_groups_file)

    # Process SMILES files for each virus
    process_virus_smiles_files(base_input_dir, virus_names, functional_groups_smarts, matched_atoms_file)

if __name__ == "__main__":
    main()
