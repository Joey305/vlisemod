import os
import csv
from pathlib import Path
import re

# List of proteins and variations to identify in PDB files
protein_variations = {
    'Protease': ['protease'],
    'Reverse Transcriptase': ['reverse transcriptase', 'transcriptase'],
    'Polymerase': ['polymerase'],
    'Helicase': ['helicase'],
    'Integrase': ['integrase'],
    'Capsid Protein': ['capsid', 'gag'],
    'Envelope Glycoprotein': ['envelope glycoprotein', 'gp120', 'env'],
    'Nucleoprotein': ['nucleoprotein', 'nucleocapsid'],
    'RNAse H': ['rnase h', 'rna h'],
    'Spike Protein': ['spike protein', 'spike'],
    'NSP Proteins': ['nsp1', 'nsp2', 'nsp3', 'nsp4', 'nsp5', 'nsp6', 'nsp7', 'nsp8', 'nsp9', 'nsp10', 'nsp11', 'nsp12', 'nsp13', 'nsp14', 'nsp15', 'nsp16'],
    'Matrix Protein': ['matrix protein'],
    'Fusion Protein': ['fusion protein'],
    'Hemagglutinin': ['hemagglutinin'],
    'Neuraminidase': ['neuraminidase'],
    'VP Proteins': ['vp1', 'vp2', 'vp3', 'vp4'],
    'Envelope Protein': ['envelope protein'],
    'Transmembrane Protein': ['transmembrane protein'],
    'Accessory Proteins': ['vpu', r'\btat\b'],  # Added regex for more precise match on TAT
    'Ribonucleoprotein': ['ribonucleoprotein', 'rnp']
}

def find_protein_in_pdb(pdb_file):
    """Check for protein strings in a PDB file and return matched proteins."""
    matched_proteins = set()
    with pdb_file.open('r') as file:
        for line in file:
            line_lower = line.lower()
            for protein, variations in protein_variations.items():
                for variant in variations:
                    if re.search(variant, line_lower):
                        matched_proteins.add(protein)
                        break  # Avoid adding the same protein multiple times
    return matched_proteins

def parse_pdb_files(virus_name, base_input_dir, output_csv):
    """Parse PDB files for a given virus and append the findings to the CSV."""
    # Define the directory where the PDB files are stored
    input_dir_path = Path(base_input_dir) / virus_name
    
    # Iterate over PDB files in the input directory
    for pdb_file in input_dir_path.glob('*.pdb'):
        pdb_id = pdb_file.stem
        
        # Find matching proteins in the PDB file
        matched_proteins = find_protein_in_pdb(pdb_file)
        
        # If proteins are found, append them to the CSV
        if matched_proteins:
            with open(output_csv, 'a', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                for protein in matched_proteins:
                    writer.writerow([virus_name, pdb_id, protein])

# Define base directories relative to the current script location
base_input_dir = './Database_DATA'
output_csv = './protein_summary.csv'

# Create output CSV file and write headers
with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['virus_name', 'pdb_id', 'protein'])

# List of viruses to process
virus_list = [
    "Human immunodeficiency virus 1",
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
    "Human papillomavirus type 6a",
    'Severe acute respiratory syndrome coronavirus 2',
]

# Run the parsing function for each virus in the list
for virus in virus_list:
    print(f"Processing virus: {virus}")
    parse_pdb_files(virus, base_input_dir, output_csv)

print("Protein extraction completed.")
