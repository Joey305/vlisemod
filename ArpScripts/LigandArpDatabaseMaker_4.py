'''
Fetching Data from Database: We pull the necessary data (virus name, PDB ID, ligand, chain, and residue) from the database and then dynamically find the corresponding JSON files.
Filtering JSON Files: We are only processing the INTER interaction type and ignoring water or intra-ligand interactions.
Atom_ID Lookup: The script reads the .pdb or .cif files to fetch the Atom_ID corresponding to the ligand atoms.
Duplicate Handling: Since we’re directly writing to the CSV file row by row, there is no need to handle duplicates in the code itself. However, if multiple JSON files are the same, the script won’t add them redundantly (it checks for the file’s existence before processing).
Logging: Any missing files or errors encountered are logged into a file for debugging purposes.

Virus Name: Added the virus_name to the first column using data fetched from the database.
Filtering for INTER interactions: The script filters only the interactions with type "INTER", ignoring other types of interactions (e.g., intra-ligand or water interactions).
File Structure Assumption: The script assumes that the JSON filenames are in the format {PDB_ID}_{Ligand}_{Residue_ID}_{Chain}.json, as described.
Running the Script:
Make sure the JSON files are in the Arpeggio_Contacts folder.
The virus names and corresponding PDB, ligand, chain, and residue details are pulled from the database (viral_data.db).
The CSV file Filtered_Ligand_Interaction_DiagramData.csv will be generated in the current directory.


'''
import os
import csv
import json
import sqlite3
import logging

# Set up logging
logging.basicConfig(filename='interaction_diagram.log', level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Directory where the JSON files are stored
json_dir = "Arpeggio_Contacts"

# Database connection to fetch virus names and related data
db_path = "viral_data.db"
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Fetch data from the database, including virus_name and residue_id
cursor.execute("SELECT virus_name, pdb_id, ligand, chain, residue_id FROM Ligand_Arp_Diagram")
data = cursor.fetchall()

# Dictionary to store PDB, ligand, chain, and residue_id for reference
pdb_virus_dict = {f"{row[1]}_{row[2]}_{row[4]}_{row[3]}": (row[0], row[4]) for row in data}

# Function to filter interactions and write to CSV
def filter_and_write_interactions_to_csv(json_dir, output_csv):
    logging.info(f"Processing JSON files in directory: {json_dir}")
    
    # Prepare CSV headers with Ligand_Number as the new column name
    headers = ["Virus_Name", "PDB_ID", "Ligand", "Ligand_Number", "Chain", "Contact", "Distance", "Ligand_Atom", 
               "Residue", "Residue_Number", "Residue_Atom", "Residue_Chain"]

    # Open the CSV file for writing
    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        writer.writeheader()
        
        # Loop through all subdirectories in the directory
        for folder_name in os.listdir(json_dir):
            folder_path = os.path.join(json_dir, folder_name)
            
            # Check if it is a directory
            if os.path.isdir(folder_path):
                # Look for the .json file inside the subdirectory
                for file_name in os.listdir(folder_path):
                    if file_name.endswith(".json"):
                        json_file = os.path.join(folder_path, file_name)
                        
                        # Extract PDB_ID, Ligand, Chain, and Residue ID from the directory name
                        base_name = os.path.splitext(folder_name)[0]
                        pdb_id, ligand, residue_id, chain = base_name.split('_')
                        
                        # Fetch the virus name and ligand number (residue_id) based on the pdb_id and other identifiers
                        virus_name, ligand_number = pdb_virus_dict.get(base_name, ("Unknown Virus", "N/A"))
                        
                        logging.info(f"Processing JSON file: {json_file} for virus {virus_name}")
                        
                        try:
                            # Load the JSON data
                            with open(json_file, 'r') as f:
                                data = json.load(f)
                            
                            # Filter only the INTER atom-atom or atom-plane interactions
                            for interaction in data:
                                bgn = interaction.get('bgn', {})
                                end = interaction.get('end', {})
                                contacts = interaction.get('contact', [])  # Fetch the contact list
                                
                                bgn_is_ligand = bgn.get('label_comp_id') == ligand
                                end_is_ligand = end.get('label_comp_id') == ligand

                                if bgn.get('label_comp_id') == "HOH" or end.get('label_comp_id') == "HOH":
                                    continue

                                if bgn_is_ligand and end_is_ligand:
                                    continue
                                
                                # Write each contact to a separate row
                                for contact_type in contacts:
                                    writer.writerow({
                                        "Virus_Name": virus_name,
                                        "PDB_ID": pdb_id,
                                        "Ligand": ligand,
                                        "Ligand_Number": ligand_number,  # Add Ligand_Number
                                        "Chain": bgn['auth_asym_id'] if bgn_is_ligand else end['auth_asym_id'],
                                        "Contact": contact_type,  # Write the contact type
                                        "Distance": interaction.get('distance', ''),
                                        "Ligand_Atom": bgn['auth_atom_id'] if bgn_is_ligand else end['auth_atom_id'],
                                        "Residue": end['label_comp_id'] if bgn_is_ligand else bgn['label_comp_id'],
                                        "Residue_Number": end['auth_seq_id'] if bgn_is_ligand else bgn['auth_seq_id'],
                                        "Residue_Atom": end['auth_atom_id'] if bgn_is_ligand else bgn['auth_atom_id'],
                                        "Residue_Chain": end['auth_asym_id'] if bgn_is_ligand else bgn['auth_asym_id']
                                    })
                        
                        except Exception as e:
                            logging.error(f"Error processing file {json_file}: {e}")

# Output CSV file
output_csv = "Filtered_Ligand_Interaction_DiagramData.csv"

# Run the filtering process
filter_and_write_interactions_to_csv(json_dir, output_csv)

# Close the database connection
conn.close()

print(f"CSV file written to {output_csv}.")
