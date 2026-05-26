import os
import requests
from bs4 import BeautifulSoup
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    filename='script_log.log',
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.DEBUG
)

# Function to fetch SMILES string from RCSB webpage
def fetch_smiles(ligand_code):
    url = f"https://www.rcsb.org/ligand/{ligand_code}"
    response = requests.get(url)
    
    if response.status_code != 200:
        logging.error(f"Failed to fetch page for ligand {ligand_code}: Status code {response.status_code}")
        return None

    soup = BeautifulSoup(response.content, 'html.parser')
    smiles_tag = soup.find('th', string='Isomeric SMILES')
    if smiles_tag:
        smiles_value = smiles_tag.find_next('td').text.strip()
        logging.info(f"Found SMILES for {ligand_code}: {smiles_value}")
        return smiles_value
    else:
        logging.warning(f"SMILES not found for {ligand_code}")
        return None

# Function to process the input .txt file and append SMILES column
def process_txt_file(input_file):
    output_file = input_file.replace('.txt', '_smiles.txt')

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = infile.readlines()
        header = lines[0].strip() + "\tSMILES\n"
        outfile.write(header)

        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) < 2:
                logging.warning(f"Skipping malformed line: {line.strip()}")
                continue
            
            ligand_code = parts[1]
            smiles = fetch_smiles(ligand_code)

            if smiles:
                outfile.write(f"{line.strip()}\t{smiles}\n")
            else:
                outfile.write(f"{line.strip()}\tN/A\n")

    logging.info(f"Processing complete. Output saved to {output_file}")

# Function to process all TXT files based on virus names
def process_virus_files(base_input_dir, virus_names):
    for virus_name in virus_names:
        txt_file = Path(base_input_dir) / 'Sorted_PDB_Files' / virus_name / 'with_ligands.txt'
        
        if not txt_file.exists():
            logging.error(f"TXT file not found for virus '{virus_name}': {txt_file}")
            continue

        logging.info(f"Processing file for virus '{virus_name}': {txt_file}")
        process_txt_file(str(txt_file))

# Main function
def main():
    base_input_dir = './Database_DATA'

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
    

    process_virus_files(base_input_dir, virus_names)

if __name__ == "__main__":
    main()
