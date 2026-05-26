import os
import csv
from pathlib import Path
import re

# List of proteins and variations to identify in CIF files
protein_variations = {
    'Protease': ['protease'],
    'Reverse Transcriptase': ['reverse transcriptase', 'transcriptase', 'TAFTIPSI'],
    'Polymerase': ['polymerase'],
    'Helicase': ['helicase'],
    'Integrase': ['integrase'],
    'Capsid Protein': ['capsid', r'\bgag\b'],
    'Envelope Protein/Glycoprotein': ['envelope glycoprotein', 'gp120', 'envelope protein', r'\benv\b', 'gp41', 'envelope', 'protein e', 'Protein E', 'Env-gp41', 'Fusion peptide', 'fusion peptide', 'AIB142', 'Aib142'],
    'Nucleoprotein': ['nucleoprotein'],
    'Nucleocapsid Protein': ['nucleocapsid'],
    'RNAse H': [r'\brnase h\b', r'\brna h\b'],
    'Spike Protein': ['spike protein', 'spike', 'spike glycoprotein' ,'Spike', 'Surface glycoprotein' 'Surface'],
    'NSP Proteins': [r'\bnsp1\b', r'\bnsp2\b', r'\bnsp3\b', r'\bnsp4\b', r'\bnsp5\b', r'\bnsp6\b', r'\bnsp7\b', r'\bnsp8\b', r'\bnsp9\b', r'\bnsp10\b', r'\bnsp11\b', r'\bnsp12\b', r'\bnsp13\b', r'\bnsp14\b', r'\bnsp15\b', r'\bnsp16\b'],
    'Matrix Protein': ['matrix protein'],
    'Fusion Protein': ['fusion protein', 'RSV fusion glycoprotein', "RSV fusion", 'rsv fusion', "Fusion glycoprotein"],
    'Hemagglutinin': ['hemagglutinin'],
    'Programmed 1 Ribosomal Frameshifting Element': ['programmed -1 ribosomal frameshifting element', 'frameshifting element', 'programmed 1 ribosomal frameshifting element'],
    'Stem-Loop RNA Element': ['Stem-Loop', 'stem loop', 'stemloop', 'Stem Loop', "5_SL4", 'Stem-Loop 4'],
    'Neuraminidase': ['neuraminidase'],
    'VP Proteins': [r'\bvp1\b', r'\bvp2\b', r'\bvp3\b', r'\bvp4\b'],
    'Transmembrane Protein': ['transmembrane protein'],
    'Accessory Proteins': [r'\bvpu\b', r'\bnef\b', r'\bNEF\b', r'\bNef\b',r'\btat\b', r'\bNef\b', r'\bvpr\b', r'\bVPR\b', r'\bVpr\b','ORF10', 'orf10', 'orf3' ,'ORF3a', 'ORF3A', 'ORF3a protein', r'\bORF\b', 'Orf6', 'ORF6', 'ORF7', 'orf7', 'ORF8', 'orf8', 'orf9', "ORF9", 'orf1', 'ORF1', 'ORF2', 'orf2', 'orf4','ORF4'],    'Ribonucleoprotein': ['ribonucleoprotein', r'\brnp\b'],
    'Exon Splicing Silencer': ['exon splicing silencer'],
    'RNA Packing Signal': ['RNA packaging signal'],
    'POL Protein': ['pol protein'],
    'Rev Protein': ['rev protein', 'Rev protein', r'\brev\b'],
    'E7 Oncoprotein' : ['oncoprotein', 'Oncoprotein',],
    "Regulatory protein E2" : ['Regulatory protein E2', r'\bE2\b'],
    'Replication protein E1' : ['Replication protein E1', r'\bE1\b'],
    'Peptide Binding Domains Complexed with MHC' : ['histocompatibility antigen ', 'histocompatibility', "HLA class II", "MHC I-peptide",],
    
    ###LEAVE COMMENTED OUT UNTIL THERE ARE NO OTHER OPTIONS FOR PAIRING AND MATCHING SPECIFIC PROTEIN STRUCTURES###
    # 'Other (Peptides & RNA Complexes)' : ['peptide', r'\bRNA\b']
    
}

def load_processed_files(output_csv):
    """Load already processed virus_name and pdb_id pairs from the output CSV."""
    processed_files = set()
    if Path(output_csv).exists():
        with open(output_csv, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            next(reader)  # Skip the header
            for row in reader:
                if len(row) >= 2:
                    virus_name, pdb_id = row[0], row[1]
                    processed_files.add((virus_name, pdb_id))
    return processed_files

def find_protein_in_cif(cif_file):
    """Check for protein strings in a CIF file and return matched proteins."""
    matched_proteins = set()
    with cif_file.open('r') as file:
        for line in file:
            line_lower = line.lower()
            for protein, variations in protein_variations.items():
                for variant in variations:
                    if re.search(variant, line_lower):
                        matched_proteins.add(protein)
                        break  # Avoid adding the same protein multiple times
    return matched_proteins

def parse_cif_files(virus_name, base_input_dir, output_csv, no_match_list, processed_files):
    """Parse CIF files for a given virus and append findings to the CSV. Log files with no matches."""
    input_dir_path = Path(base_input_dir) / virus_name
    
    # Iterate over CIF files in the input directory
    for cif_file in input_dir_path.glob('*.cif'):
        pdb_id = cif_file.stem

        # Skip if the file has already been processed
        if (virus_name, pdb_id) in processed_files:
            continue
        
        # Find matching proteins in the CIF file
        matched_proteins = find_protein_in_cif(cif_file)
        
        if matched_proteins:
            # If proteins are found, append them to the CSV
            with open(output_csv, 'a', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                for protein in matched_proteins:
                    writer.writerow([virus_name, pdb_id, protein])
        else:
            # If no match, append the file to the no match list
            no_match_list.append(cif_file.name)

# Define base directories relative to the current script location
base_input_dir = './Arpeggio_Data'
output_csv = './protein_summaryCIF.csv'
no_match_file = './no_protein_match_files.txt'

# Create or load the output CSV file and load already processed virus_name and pdb_id pairs
processed_files = load_processed_files(output_csv)

# List to track files with no matches
no_match_list = []

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
    parse_cif_files(virus, base_input_dir, output_csv, no_match_list, processed_files)

# Write the list of files with no protein match to a text file
with open(no_match_file, 'w', encoding='utf-8') as f:
    for no_match in no_match_list:
        f.write(f"{no_match}\n")

print("Protein extraction completed.")
print(f"Files with no matches written to {no_match_file}.")
