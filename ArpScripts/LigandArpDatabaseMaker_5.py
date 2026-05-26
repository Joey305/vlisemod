import os
import csv
import logging

# Set up logging
logging.basicConfig(filename='atom_id_log.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# CSV input and output file paths
csv_file = "Filtered_Ligand_Interaction_DiagramData.csv"
output_csv_file = "Updated_Ligand_Interaction_DiagramData.csv"

# Function to extract Atom_Number (Atom_ID) from a .pdb file based on CSV row information
def extract_atom_number_from_pdb(pdb_file, ligand_atom, ligand_code, chain, ligand_number):
    atom_number = None
    logging.info(f"Extracting Atom_ID from {pdb_file} for Ligand_Atom: {ligand_atom}, Ligand_Code: {ligand_code}, Chain: {chain}, Ligand_Number: {ligand_number}")

    try:
        with open(pdb_file, 'r') as pdb:
            for line in pdb:
                if line.startswith('HETATM'):  # Only parse HETATM lines in PDB file
                    pdb_atom_number = line[6:11].strip()  # Atom serial number
                    pdb_atom_name = line[12:16].strip()  # Atom name
                    pdb_ligand_code = line[17:20].strip()  # Ligand code
                    pdb_chain = line[21].strip()  # Chain ID
                    pdb_residue_number = line[22:26].strip()  # Residue number
                    
                    # Log what is being processed from the PDB file
                    logging.info(f"Processing PDB line: Atom_Name: {pdb_atom_name}, Ligand_Code: {pdb_ligand_code}, Chain: {pdb_chain}, Residue_Number: {pdb_residue_number}, Atom_Number: {pdb_atom_number}")

                    # Match the ligand atom, ligand code, chain, and residue number (which corresponds to Ligand_Number in CSV)
                    if pdb_atom_name == ligand_atom and pdb_ligand_code == ligand_code and pdb_chain == chain and pdb_residue_number == ligand_number:
                        logging.info(f"Match found! Atom_Number: {pdb_atom_number}")
                        atom_number = pdb_atom_number
                        return atom_number  # Return the matched atom number
    
    except Exception as e:
        logging.error(f"Error reading {pdb_file}: {e}")

    logging.warning(f"Atom_ID not found for Ligand_Atom: {ligand_atom}, Ligand_Code: {ligand_code}, Chain: {chain}, Ligand_Number: {ligand_number} in {pdb_file}")
    return atom_number  # Return None if not found

# Function to process each row and append to the output CSV
def process_csv_row_by_row(input_csv, output_csv, pdb_base_dir):
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        
        # Get the fieldnames and insert 'Atom_ID' right after 'Ligand_Atom'
        fieldnames = reader.fieldnames[:]
        atom_id_index = fieldnames.index('Ligand_Atom') + 1  # Find the position after 'Ligand_Atom'
        fieldnames.insert(atom_id_index, 'Atom_ID')  # Insert 'Atom_ID' after 'Ligand_Atom'
        
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        # Write the header to the output CSV
        writer.writeheader()

        for row in reader:
            pdb_id = row['PDB_ID']
            virus_name = row['Virus_Name']
            ligand_atom = row['Ligand_Atom']
            ligand_code = row['Ligand']
            chain = row['Chain']
            ligand_number = str(row['Ligand_Number'])  # Ligand_Number corresponds to Residue_Number
            
            # Construct the path to the PDB file based on the virus name and PDB ID
            pdb_file = os.path.join(pdb_base_dir, virus_name, f"{pdb_id}.pdb")

            # Check if the PDB file exists
            if not os.path.exists(pdb_file):
                logging.warning(f'PDB file not found for PDB_ID: {pdb_id} in Virus directory: {virus_name}')
                continue  # Skip this row if the PDB file is not found

            # Extract Atom_Number (Atom_ID) from the PDB file
            atom_id = extract_atom_number_from_pdb(pdb_file, ligand_atom, ligand_code, chain, ligand_number)
            
            if atom_id:  # Only write the row if Atom_ID is found
                row['Atom_ID'] = atom_id
                writer.writerow(row)
            else:
                logging.warning(f'Skipping row due to missing Atom_ID for PDB_ID: {pdb_id}, Ligand_Atom: {ligand_atom}, Ligand_Code: {ligand_code}, Chain: {chain}, Ligand_Number: {ligand_number}')

# Directory containing the virus-specific PDB files
pdb_base_directory = "Database_Data"  # This is the base directory containing virus directories

# Process the CSV file row by row and write to a new CSV file
process_csv_row_by_row(csv_file, output_csv_file, pdb_base_directory)

print(f"CSV file has been processed row-by-row. Output written to {output_csv_file}. Check 'atom_id_log.log' for warnings/errors.")
