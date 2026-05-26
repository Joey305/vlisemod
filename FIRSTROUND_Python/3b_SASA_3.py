"""
# PyMOL Automation Script for Processing Viral PDB Data

## Overview
This script automates the processing and visualization of viral protein-ligand interactions using PyMOL. It identifies atoms within specified distances to water molecules and protein atoms, creates a CSV output for further analysis, and generates PyMOL sessions where these atoms are visualized based on their properties (e.g., solvent exposure, proximity to binding pocket). 

The script focuses on:
1. Extracting ligand atoms from the PDB files.
2. Filtering ligand atoms by proximity to water molecules and proteins.
3. Identifying solvent-exposed and distal atoms.
4. Generating binding pocket atoms.
5. Automatically creating a PyMOL session to visualize the ligand, protein, and water interactions, colored by atom type.

## Required Setup
0. STRIP THE TIMESTAMP OFF THE CSV FILE GENERATED DURING STEP 2
### 1. Software Requirements
- **Python 3.x**: Ensure that Python is installed on your system.
- **PyMOL**: PyMOL must be installed and the path to the executable set correctly in the script. The default path to the PyMOL executable is set as:

  pymol_executable = r"C:\ProgramData\pymol\PyMOLWin.exe"

  Update this path according to your system configuration if it differs.

2. Directory and File Structure
The script expects a specific directory structure to function properly:

Database_DATA/: This folder should contain a subdirectory for each virus type, named after the virus (e.g., Human immunodeficiency virus 1). Each of these subdirectories should contain PDB files with structures of the virus.
Database_DATA/Sorted_PDB_Files/: Inside this directory, there should be a CSV file for each virus, containing the PDB IDs and corresponding ligand codes for each structure. The CSV file should be named using the virus name (e.g., with_ligands_test.csv).
3. Input CSV File
Each CSV file should have the following structure:

Copy code
PDB ID,Ligand Codes
1A8G,2Z4
This structure helps the script identify which PDB files to process and which ligands to extract for each PDB file.

4. Output Structure
The script generates outputs in the Processed_PDBs/ directory. For each PDB-ligand combination, a subdirectory is created where various CSV files are saved, along with the generated PyMOL session.

5. Running the Script
To run the script, navigate to the directory containing the script and use the following command:

bash
Copy code
python 3_SASA_3.py
The main workflow function main_workflow() will automatically process all the PDB files listed in the input CSV, generate the required CSVs, and create the PyMOL visualization.

Detailed Script Explanation
Main Functions
setup_output_directory()

Creates the necessary output directory for each virus and PDB-ligand pair.
extract_ligand_atoms()

Extracts ligand atoms from the PDB file and writes them to a CSV file.
process_pdb_file_by_water()

Filters ligand atoms by proximity to water molecules, saving the result to a CSV.
filter_by_water_vs_protein()

Compares ligand atoms' proximity to water molecules and protein atoms to identify solvent-exposed atoms.
create_distal_atoms_csv()

Identifies distal ligand atoms that have no nearby protein atoms.
identify_binding_pocket()

Identifies binding pocket atoms that are near the ligand and saves them to a CSV.
write_pymol_script()

Generates a PyMOL script to visualize the ligand, protein, and water molecules. Ligand atoms are shown as sticks and colored by atom type (e.g., carbon, oxygen, nitrogen, sulfur). Other key atoms (distal, solvent-exposed) are also highlighted using spheres and appropriate color schemes.
execute_pymol_script()

Executes the generated PyMOL script, creating a saved PyMOL session for further analysis.
Output CSV Files
_ligand.csv: Contains the extracted ligand atoms.
_ligand_water_distances.csv: Lists the distances between ligand atoms and nearby water molecules.
_solvent_exposed_atoms.csv: Contains solvent-exposed ligand atoms.
_distalAtoms.csv: Lists distal ligand atoms with no nearby protein atoms.
_receptor_binding_pocket.csv: Contains atoms within the receptor's binding pocket near the ligand.
PyMOL Visualization
The generated PyMOL session contains:

Ligand_Object: The ligand displayed as sticks and colored by atom type.
Filtered_Water_Protein: Solvent-exposed atoms near water molecules.
Distal_Atoms: Distal ligand atoms shown as spheres.
Binding_Pocket: Binding pocket atoms highlighted as spheres.
The protein is displayed in cartoon format, and various atoms of interest are highlighted with different color schemes.

Notes:
Ensure that PyMOL is properly installed and callable from the command line. If PyMOL cannot be found, update the path to the PyMOL executable in the script.
This script assumes the PDB files contain valid ligand and water atom entries. Ensure the PDB files are correct and complete for proper functionality. """
  

  

  

import os
import csv
import shutil
import subprocess
from pathlib import Path
from math import sqrt

# Path to the PyMOL executable
pymol_executable = r"C:\ProgramData\pymol\PyMOLWin.exe"

# # Virus name variable at the top of the script
# virus_name = 'Severe acute respiratory syndrome coronavirus 2'

# Set the number of nearby protein atoms threshold
protein_atom_threshold = 2

# Helper function to set up the output directory
def setup_output_directory(base_output_dir, virus_name, pdb_id, ligand):
    virus_dir = Path(base_output_dir) / virus_name
    virus_dir.mkdir(parents=True, exist_ok=True)

    folder_name = f"{pdb_id}_{ligand}"
    output_dir = virus_dir / folder_name
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Created directory: {output_dir}")
    return output_dir

# Helper function to log information to the log.txt file
def log_info(log_file, message):
    with open(log_file, 'a') as log:
        log.write(message + '\n')

# Helper function to extract all ligand atoms and save them in a CSV
def extract_ligand_atoms(pdb_file, ligand, output_dir):
    ligand_atoms = []
    log_file = os.path.join(output_dir, "log.txt")

    with open(pdb_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("HETATM"):
                residue_name = line[17:20].strip()
                if residue_name == ligand:
                    atom_id = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    chain = line[21].strip()
                    sequence_id = line[22:26].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atom_type = line[76:78].strip()

                    atom_entry = {
                        "PDB File": pdb_file,
                        "Ligand": residue_name,
                        "Chain": chain,
                        "Sequence ID": atom_id,
                        "Exact Atom": atom_name,
                        "Atom Type": atom_type,
                        "X": x,
                        "Y": y,
                        "Z": z
                    }
                    ligand_atoms.append(atom_entry)

                    # Log the extracted atom details
                    log_info(log_file, f"Extracted ligand atom: {atom_name} (ID: {atom_id}) at coordinates ({x}, {y}, {z})")

    output_csv_ligand_atoms = os.path.join(output_dir, f"{Path(pdb_file).stem}_{ligand}_ligand.csv")
    with open(output_csv_ligand_atoms, 'w', newline='') as csvfile:
        fieldnames = ["PDB File", "Ligand", "Chain", "Sequence ID", "Exact Atom", "Atom Type", "X", "Y", "Z"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ligand_atoms)

    print(f"CSV (ligand atoms) created: {output_csv_ligand_atoms}")
    return output_csv_ligand_atoms

# Helper function to process each PDB file and filter by water proximity
def process_pdb_file_by_water(pdb_file, ligand, output_dir, water_distance_threshold=5.0):
    atom_data = []
    ligand_atoms = []
    water_atoms = []
    log_file = os.path.join(output_dir, "log.txt")

    with open(pdb_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("HETATM"):
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain = line[21].strip()
                sequence_id = line[22:26].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom_type = line[76:78].strip()

                atom_entry = {
                    "PDB File": pdb_file,
                    "Ligand": residue_name,
                    "Chain": chain,
                    "Sequence ID": atom_id,
                    "Exact Atom": atom_name,
                    "Atom Type": atom_type,
                    "X": x,
                    "Y": y,
                    "Z": z
                }

                if residue_name == ligand:
                    ligand_atoms.append(atom_entry)
                elif residue_name == "HOH":
                    water_atoms.append(atom_entry)

    for ligand_atom in ligand_atoms:
        for water_atom in water_atoms:
            dx = (ligand_atom["X"] - water_atom["X"]) ** 2
            dy = (ligand_atom["Y"] - water_atom["Y"]) ** 2
            dz = (ligand_atom["Z"] - water_atom["Z"]) ** 2
            distance = sqrt(dx + dy + dz)

            if distance <= water_distance_threshold:
                atom_data.append({
                    "PDB File": ligand_atom["PDB File"],
                    "Ligand": ligand_atom["Ligand"],
                    "Chain": ligand_atom["Chain"],
                    "Sequence ID": ligand_atom["Sequence ID"],
                    "Exact Atom": ligand_atom["Exact Atom"],
                    "Atom Type": ligand_atom["Atom Type"],
                    "X": ligand_atom["X"],
                    "Y": ligand_atom["Y"],
                    "Z": ligand_atom["Z"],
                    "Water Chain": water_atom["Chain"],
                    "Water Sequence ID": water_atom["Sequence ID"],
                    "Water X": water_atom["X"],
                    "Water Y": water_atom["Y"],
                    "Water Z": water_atom["Z"],
                    "Distance": distance
                })

                # Log the ligand atom that is near a water molecule
                log_info(log_file, f"Ligand atom {ligand_atom['Exact Atom']} (ID: {ligand_atom['Sequence ID']}) is {distance:.2f} Å from water molecule at coordinates ({water_atom['X']}, {water_atom['Y']}, {water_atom['Z']})")

    output_csv_by_water = os.path.join(output_dir, f"{Path(pdb_file).stem}_{ligand}_ligand_water_distances.csv")
    with open(output_csv_by_water, 'w', newline='') as csvfile:
        fieldnames = [
            "PDB File", "Ligand", "Chain", "Sequence ID", "Exact Atom", "Atom Type", "X", "Y", "Z",
            "Water Chain", "Water Sequence ID", "Water X", "Water Y", "Water Z", "Distance"
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(atom_data)

    print(f"CSV (by water) created: {output_csv_by_water}")
    return output_csv_by_water

# Function to further filter atoms by proximity to protein and retain unique atoms
def filter_by_water_vs_protein(csv_file, pdb_file, output_dir, ligand, protein_distance_threshold=4.15):
    unique_atoms = {}
    protein_atoms = []
    log_file = os.path.join(output_dir, "log.txt")

    # Parse the CSV file (ligand_water_distances.csv) to extract water-proximity data
    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ligand_atom_key = (row['Chain'], row['Sequence ID'], row['Exact Atom'])
            if ligand_atom_key not in unique_atoms:
                unique_atoms[ligand_atom_key] = {
                    "PDB File": row["PDB File"],
                    "Ligand": row["Ligand"],
                    "Chain": row["Chain"],
                    "Sequence ID": row["Sequence ID"],
                    "Exact Atom": row["Exact Atom"],
                    "Atom Type": row["Atom Type"],
                    "X": row["X"],
                    "Y": row["Y"],
                    "Z": row["Z"],
                    "Water Proximity Count": 0,
                    "Protein Proximity Count": 0
                }
            unique_atoms[ligand_atom_key]["Water Proximity Count"] += 1

    # Parse the PDB file for protein atoms
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                atom_id = int(line[6:11].strip())
                chain = line[21].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                protein_atoms.append({
                    "Atom ID": atom_id,
                    "Chain": chain,
                    "X": x,
                    "Y": y,
                    "Z": z
                })

    # Compare each ligand atom to protein atoms and calculate proximity
    for atom_key, atom_data in unique_atoms.items():
        ligand_x = float(atom_data["X"])
        ligand_y = float(atom_data["Y"])
        ligand_z = float(atom_data["Z"])

        for protein_atom in protein_atoms:
            dx = (ligand_x - protein_atom["X"]) ** 2
            dy = (ligand_y - protein_atom["Y"]) ** 2
            dz = (ligand_z - protein_atom["Z"]) ** 2
            distance = sqrt(dx + dy + dz)

            if distance <= protein_distance_threshold:
                atom_data["Protein Proximity Count"] += 1

        # Log the proximity counts for water and protein for each atom
        log_info(log_file, f"Ligand atom {atom_data['Exact Atom']} (ID: {atom_data['Sequence ID']}) has {atom_data['Water Proximity Count']} water molecules within 5 Å and {atom_data['Protein Proximity Count']} protein atoms within 4.15 Å.")

    # Retain only atoms that have more nearby water molecules than protein atoms
    solvent_exposed_atoms = [
        {
            "PDB File": atom_data["PDB File"],
            "Ligand": atom_data["Ligand"],
            "Chain": atom_data["Chain"],
            "Sequence ID": atom_data["Sequence ID"],
            "Exact Atom": atom_data["Exact Atom"],
            "Atom Type": atom_data["Atom Type"],
            "X": atom_data["X"],
            "Y": atom_data["Y"],
            "Z": atom_data["Z"]
        }
        for atom_data in unique_atoms.values()
        if atom_data["Water Proximity Count"] > atom_data["Protein Proximity Count"]
    ]

    # Write the solvent-exposed atoms to a CSV
    output_csv_by_protein = os.path.join(output_dir, f"{Path(pdb_file).stem}_{ligand}_solvent_exposed_atoms.csv")
    with open(output_csv_by_protein, 'w', newline='') as csvfile:
        fieldnames = ["PDB File", "Ligand", "Chain", "Sequence ID", "Exact Atom", "Atom Type", "X", "Y", "Z"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(solvent_exposed_atoms)

    print(f"CSV (solvent-exposed atoms) created: {output_csv_by_protein}")
    return output_csv_by_protein

# Function to create distal atoms CSV from ligand atoms
def create_distal_atoms_csv(ligand_csv, pdb_file, output_dir, ligand, protein_distance_threshold=4.15):
    distal_atoms = []
    ligand_atoms = []
    protein_atoms = []
    log_file = os.path.join(output_dir, "log.txt")

    with open(ligand_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ligand_atoms.append(row)

    with open(pdb_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                atom_id = int(line[6:11].strip())
                chain = line[21].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                protein_atoms.append({
                    "Atom ID": atom_id,
                    "Chain": chain,
                    "X": x,
                    "Y": y,
                    "Z": z
                })

    for ligand_atom in ligand_atoms:
        nearby_protein_atoms = 0

        ligand_x = float(ligand_atom["X"])
        ligand_y = float(ligand_atom["Y"])
        ligand_z = float(ligand_atom["Z"])

        for protein_atom in protein_atoms:
            dx = (ligand_x - protein_atom["X"]) ** 2
            dy = (ligand_y - protein_atom["Y"]) ** 2
            dz = (ligand_z - protein_atom["Z"]) ** 2
            distance = sqrt(dx + dy + dz)

            if distance <= protein_distance_threshold:
                nearby_protein_atoms += 1

        if nearby_protein_atoms == 0:
            distal_atoms.append(ligand_atom)

            # Log the distal atom details
            log_info(log_file, f"Ligand atom {ligand_atom['Exact Atom']} (ID: {ligand_atom['Sequence ID']}) retained as distal atom because no protein atoms are within 4.15 Å.")

    output_csv_distal_atoms = os.path.join(output_dir, f"{Path(pdb_file).stem}_{ligand}_distalAtoms.csv")
    with open(output_csv_distal_atoms, 'w', newline='') as csvfile:
        fieldnames = list(ligand_atoms[0].keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(distal_atoms)

    print(f"CSV (distal atoms) created: {output_csv_distal_atoms}")
    return output_csv_distal_atoms


# Helper function to identify atoms within the binding pocket and save them
def identify_binding_pocket(pdb_file, ligand_csv, output_dir, binding_distance_threshold=5.0):
    binding_pocket_atoms = []
    ligand_atoms = []

    # Load ligand atoms from the ligand CSV
    with open(ligand_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ligand_atoms.append({
                "X": float(row["X"]),
                "Y": float(row["Y"]),
                "Z": float(row["Z"])
            })

    # Process PDB file to identify atoms near the ligand atoms
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain = line[21].strip()
                residue_sequence_id = int(line[22:26].strip())  # Sequence ID of the residue
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom_type = line[76:78].strip()

                # Compare distance to each ligand atom
                for ligand_atom in ligand_atoms:
                    dx = (x - ligand_atom["X"]) ** 2
                    dy = (y - ligand_atom["Y"]) ** 2
                    dz = (z - ligand_atom["Z"]) ** 2
                    distance = sqrt(dx + dy + dz)

                    if distance <= binding_distance_threshold:
                        binding_pocket_atoms.append({
                            "PDB File": pdb_file,
                            "Residue": f"{residue_name}{residue_sequence_id}",
                            "Chain": chain,
                            "Atom ID": atom_id,  # Ensure we use Atom ID here
                            "Exact Atom": atom_name,
                            "Atom Type": atom_type,
                            "X": x,
                            "Y": y,
                            "Z": z
                        })
                        break  # Stop once we find a matching atom

    # Write the binding pocket atoms to a CSV
    output_csv_binding_pocket = os.path.join(output_dir, f"{Path(pdb_file).stem}_receptor_binding_pocket.csv")
    with open(output_csv_binding_pocket, 'w', newline='') as csvfile:
        fieldnames = ["PDB File", "Residue", "Chain", "Atom ID", "Exact Atom", "Atom Type", "X", "Y", "Z"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(binding_pocket_atoms)

    print(f"CSV (binding pocket atoms) created: {output_csv_binding_pocket}")
    return output_csv_binding_pocket


# Function to write the PyMOL script
def write_pymol_script(pdb_file, hydrated_csv_file, distal_csv_file, solvent_exposed_csv_file, binding_pocket_csv_file, ligand_code, output_script):
    hydrated_atoms = set()
    distal_atoms = set()
    solvent_exposed_atoms = set()
    binding_pocket_atoms = set()

    # Load hydrated atoms
    with open(hydrated_csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            atom_id = row['Sequence ID']
            chain = row['Chain']
            hydrated_atoms.add((atom_id, chain))

    # Load distal atoms
    with open(distal_csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            atom_id = row['Sequence ID']
            chain = row['Chain']
            distal_atoms.add((atom_id, chain))

    # Load solvent-exposed atoms
    with open(solvent_exposed_csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            atom_id = row['Sequence ID']
            chain = row['Chain']
            solvent_exposed_atoms.add((atom_id, chain))

    # Load binding pocket atoms
    with open(binding_pocket_csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            atom_id = row['Atom ID']
            chain = row['Chain']
            binding_pocket_atoms.add((atom_id, chain))

    with open(output_script, 'w') as script:
        # Load the PDB file
        script.write(f"load {pdb_file}\n")

        # Create selections for hydrated atoms
        script.write("create Hydrated_Atoms, none\n")
        for atom_id, chain in hydrated_atoms:
            selection = f"chain {chain} and id {atom_id}"
            script.write(f"select temp, {selection}\n")
            script.write(f"create Hydrated_Atoms, Hydrated_Atoms or temp\n")
            script.write(f"delete temp\n")

        # Create selections for distal atoms
        script.write("create Distal_Atoms, none\n")
        for atom_id, chain in distal_atoms:
            selection = f"chain {chain} and id {atom_id}"
            script.write(f"select temp, {selection}\n")
            script.write(f"create Distal_Atoms, Distal_Atoms or temp\n")
            script.write(f"delete temp\n")

        # Create selections for solvent-exposed atoms
        script.write("create Filtered_Water_Protein, none\n")
        for atom_id, chain in solvent_exposed_atoms:
            selection = f"chain {chain} and id {atom_id}"
            script.write(f"select temp, {selection}\n")
            script.write(f"create Filtered_Water_Protein, Filtered_Water_Protein or temp\n")
            script.write(f"delete temp\n")

        # Create selections for binding pocket atoms
        script.write("create Binding_Pocket, none\n")
        for atom_id, chain in binding_pocket_atoms:
            selection = f"chain {chain} and id {atom_id}"
            script.write(f"select temp, {selection}\n")
            script.write(f"create Binding_Pocket, Binding_Pocket or temp\n")
            script.write(f"delete temp\n")

        # Create a new object for the ligand
        script.write(f"select {ligand_code}_Ligand, resn {ligand_code}\n")
        script.write(f"create Ligand_Object, {ligand_code}_Ligand\n")  # Create a new object for the ligand
        script.write(f"show sticks, Ligand_Object\n")  # Show the ligand as sticks
        # Color by atom type
        script.write(f"color gray, Ligand_Object and elem C\n")  # Color carbon atoms gray
        script.write(f"color red, Ligand_Object and elem O\n")     # Color oxygen atoms red
        script.write(f"color blue, Ligand_Object and elem N\n")    # Color nitrogen atoms blue
        script.write(f"color yellow, Ligand_Object and elem S\n")   # Color sulfur atoms yellow

        # Show spheres for Distal_Atoms
        script.write("show spheres, Distal_Atoms\n")
        script.write("color blue, Distal_Atoms\n")

        # Show spheres for Filtered_Water_Protein
        script.write("show spheres, Filtered_Water_Protein\n")
        script.write("color firebrick, Filtered_Water_Protein\n")

        # Show sticks for Hydrated_Atoms
        script.write("show sticks, Hydrated_Atoms\n")
        script.write("color cyan, Hydrated_Atoms\n")

        # Show cartoon for the protein structure and hide lines
        script.write("show cartoon, all\n")
        script.write("hide lines, all\n")

        # Show spheres for Binding_Pocket atoms
        script.write("show spheres, Binding_Pocket\n")
        script.write("color yellow, Binding_Pocket\n")

        # Save the PyMOL session
        session_file = Path(output_script).with_suffix('.pse')  # Full path for session file
        script.write(f"save {session_file}\n")



def execute_pymol_script(pymol_script_path):
    try:
        subprocess.run([pymol_executable, "-cq", pymol_script_path], check=True)
        print(f"PyMOL session created and saved using {pymol_script_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing PyMOL script: {e}")

# Function to process multiple viruses
def main_workflow(virus_names, base_output_dir, water_distance_threshold=5.0, protein_distance_threshold=4.15, protein_atom_threshold=2, binding_distance_threshold=5.0):
    script_dir = Path(__file__).parent.resolve()

    for virus_name in virus_names:
        print(f"Processing virus: {virus_name}")
        
        # Dynamically locate the latest .txt file starting with "with_" in the output directory for the current virus
        sorted_pdb_dir = script_dir / f"Database_DATA/Sorted_PDB_Files/{virus_name}"
        txt_files = list(sorted_pdb_dir.glob("with_*.txt"))
        
        if not txt_files:
            print(f"Error: No .txt files starting with 'with_' found in {sorted_pdb_dir}. Skipping virus: {virus_name}")
            continue

        # Select the latest .txt file based on the timestamp in the filename
        latest_txt_file = max(txt_files, key=os.path.getctime)
        print(f"Using input file: {latest_txt_file} for virus: {virus_name}")

        # Open the latest .txt file
        with open(latest_txt_file, 'r') as file:
            lines = file.readlines()

        for line in lines:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue

            # Split the line by tabs
            parts = line.split('\t')
            if len(parts) < 2:
                print(f"Skipping malformed line: {line}")
                continue
            
            pdb_id = parts[0]
            ligand_codes = parts[1]
            ligands = ligand_codes.split(', ')  # Split ligands by comma and space if there are multiple ligands

            pdb_file_path = script_dir / f"Database_DATA/{virus_name}/{pdb_id}.pdb"

            if not pdb_file_path.exists():
                print(f"Error: {pdb_file_path} not found for virus: {virus_name}. Skipping.")
                continue

            # Process each ligand in the ligands list
            for ligand in ligands:
                output_dir = setup_output_directory(base_output_dir, virus_name, pdb_id, ligand)

                shutil.copy(pdb_file_path, output_dir)
                print(f"Copied PDB file: {pdb_file_path} to {output_dir}")

                # Generate the four CSV files
                csv_output_ligand_atoms = extract_ligand_atoms(pdb_file_path, ligand, output_dir)
                csv_output_by_water = process_pdb_file_by_water(pdb_file_path, ligand, output_dir, water_distance_threshold)
                csv_output_by_protein = filter_by_water_vs_protein(csv_output_by_water, pdb_file_path, output_dir, ligand, protein_distance_threshold)
                csv_output_distal_atoms = create_distal_atoms_csv(csv_output_ligand_atoms, pdb_file_path, output_dir, ligand, protein_distance_threshold)

                # Identify the binding pocket and generate the binding pocket CSV
                csv_output_binding_pocket = identify_binding_pocket(pdb_file_path, csv_output_ligand_atoms, output_dir, binding_distance_threshold)

                # Create and execute the PyMOL script
                pymol_script_path = os.path.join(output_dir, f"{Path(pdb_file_path).stem}_pymol_session.pml")
                write_pymol_script(pdb_file_path, csv_output_by_water, csv_output_distal_atoms, csv_output_by_protein, csv_output_binding_pocket, ligand, pymol_script_path)

                execute_pymol_script(pymol_script_path)

# Set output directory and list of viruses
base_output_dir = "./Processed_PDBs/"
virus_names = [
    # "Severe acute respiratory syndrome coronavirus 2",
    # "Human immunodeficiency virus 1"
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
    # "Human immunodeficiency virus 1",


]

# Run the workflow for multiple viruses
main_workflow(virus_names, base_output_dir, water_distance_threshold=5.0, protein_distance_threshold=4.15, protein_atom_threshold=2, binding_distance_threshold=5.0)