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

# Function to identify functional groups in a molecule
def identify_functional_groups(smiles, functional_groups_smarts):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.error(f"Failed to parse SMILES: {smiles}")
        return []

    functional_groups = []
    for group_name, smarts in functional_groups_smarts.items():
        fg_mol = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(fg_mol):
            functional_groups.append(group_name)
    
    return functional_groups

# Function to process the SMILES .txt file
def process_smiles_file(input_file, functional_groups_smarts):
    output_file = input_file.replace('.txt', '_functional_groups_GROUPED.txt')

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        # Read and write header with additional Functional Groups column
        header = next(reader)
        header.append('Functional Groups')
        writer.writerow(header)

        for row in reader:
            smiles = row[-1]  # Assuming SMILES is in the last column
            functional_groups = identify_functional_groups(smiles, functional_groups_smarts)
            row.append(", ".join(functional_groups) if functional_groups else "None")
            writer.writerow(row)
            logging.info(f"Processed ligand with SMILES: {smiles}")

    logging.info(f"Processing complete. Output saved to {output_file}")

# Function to process all SMILES files based on virus names
def process_virus_smiles_files(base_input_dir, virus_names, functional_groups_smarts):
    for virus_name in virus_names:
        smiles_file = Path(base_input_dir) / 'Sorted_PDB_Files' / virus_name / 'with_ligands_smiles.txt'
        
        if not smiles_file.exists():
            logging.error(f"SMILES file not found for virus '{virus_name}': {smiles_file}")
            continue

        logging.info(f"Processing SMILES file for virus '{virus_name}': {smiles_file}")
        process_smiles_file(str(smiles_file), functional_groups_smarts)

# Main function
def main():
    base_input_dir = './Database_DATA'
    functional_groups_file = 'functional_groups.txt'  # Path to your functional groups file

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

    process_virus_smiles_files(base_input_dir, virus_names, functional_groups_smarts)

if __name__ == "__main__":
    main()
