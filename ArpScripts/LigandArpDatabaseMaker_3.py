import sqlite3
import subprocess
import os

# Step 1: Connect to your viral_data.db database and fetch relevant data
db_path = "viral_data.db"  # Replace with the path to your database
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Fetch data from the Ligand_Arp_Diagram.csv table (or equivalent)
cursor.execute("SELECT virus_name, pdb_id, ligand, chain, residue_id FROM Ligand_Arp_Diagram")
data = cursor.fetchall()

# Step 2: Define a function to run pdbe-arpeggio
def run_arpeggio(pdb_file, ligand, chain, residue_id, output_dir="Arpeggio_Contacts"):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Name the output file based on the naming convention
    output_file = os.path.join(output_dir, f"{os.path.basename(pdb_file).split('.')[0]}_{ligand}_{residue_id}_{chain}.json")

    # Quote the file paths to handle spaces in directory names
    pdb_file = f'"{pdb_file}"'  # Keep quotes for command usage
    output_file_path = os.path.join(output_dir, f"{os.path.basename(pdb_file).split('.')[0]}_{ligand}_{residue_id}_{chain}.json")

    # Specify the ligand selection (e.g., /A/508/)
    selection = f"/{chain}/{residue_id}/"

    # Run the pdbe-arpeggio command with the selection
    command = f'pdbe-arpeggio -s {selection} -o "{output_file_path}" {pdb_file}'
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Arpeggio analysis for {pdb_file}, ligand {ligand}, chain {chain}, residue {residue_id} completed.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to run Arpeggio for {pdb_file} and ligand {ligand}: {e}")
    except OSError as e:
        print(f"OS error occurred: {e}")

# Step 3: Iterate through each row from the database query and run Arpeggio
for row in data:
    virus_name, pdb_id, ligand, chain, residue_id = row

    # Check for either .pdb or .cif file, prioritize CIF
    pdb_file = f"Arpeggio_DATA/{virus_name}/{pdb_id}.cif" if os.path.exists(f"Arpeggio_DATA/{virus_name}/{pdb_id}.cif") else f"Arpeggio_DATA/{virus_name}/{pdb_id}.pdb"
    
    if os.path.exists(pdb_file):
        print(f"Processing {pdb_id} with ligand {ligand}, chain {chain}, residue {residue_id}")
        run_arpeggio(pdb_file, ligand, chain, residue_id)
    else:
        print(f"File {pdb_file} not found!")

# Step 4: Close the database connection when done
conn.close()
