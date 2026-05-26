import os
import csv
import subprocess
from pathlib import Path
from rdkit import Chem

# Path to the PyMOL executable (Modify this to your correct PyMOL path)
pymol_executable = r"C:\ProgramData\pymol\PyMOLWin.exe"

# Helper function to set up the output directory
def setup_output_directory(base_output_dir, virus_name, pdb_id, ligand, chain):
    output_dir = Path(base_output_dir) / virus_name / f"{pdb_id}_{ligand}_{chain}"
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Created directory: {output_dir}")
    return output_dir

# Function to execute the PyMOL script
def execute_pymol_script(pymol_script_path):
    try:
        subprocess.run([pymol_executable, "-cq", pymol_script_path], check=True)
        print(f"PyMOL session created and saved using {pymol_script_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing PyMOL script: {e}")

# Main function to process each entry
def process_pdb_ligand_chain(pdb_id, ligand, chain, pdb_file, functional_group_info, output_dir):
    functional_groups = {}
    
    # Read the functional groups for the specific PDB, ligand, and chain
    for fg_name, atom_id, exact_atom, atom_type in functional_group_info:
        if fg_name not in functional_groups:
            functional_groups[fg_name] = []
        functional_groups[fg_name].append((atom_id, exact_atom, atom_type))
    
    # Set up the output directory for this combo
    output_script = output_dir / f"{pdb_id}_{ligand}_{chain}_pymol_session.pml"
    
    # Write and execute the PyMOL script
    write_pymol_script(pdb_file, output_script, ligand, chain, functional_groups)
    execute_pymol_script(output_script)

    
    # Write and execute the PyMOL script
    write_pymol_script(pdb_file, output_script, ligand, chain, functional_groups)
    execute_pymol_script(output_script)

# Function to process the functional groups file from Sorted_PDB_Files directory
def process_functional_groups_file(virus_name, base_output_dir):
    sorted_pdb_dir = Path(f"./Database_DATA/Sorted_PDB_Files/{virus_name}")
    functional_groups_file = sorted_pdb_dir / "with_ligands_smiles_functional_groups-reorganized.txt"

    if not functional_groups_file.exists():
        print(f"Functional groups file not found for virus: {virus_name}")
        return

    with open(functional_groups_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        next(reader)  # Skip header

        current_pdb = None
        current_ligand = None
        current_chain = None
        pdb_file = None
        functional_group_info = []

        for row in reader:
            # Unpack all seven columns
            pdb_id, ligand, chain, functional_group, atom_id, exact_atom, atom_type = row[:7]
            atom_id = int(atom_id)

            pdb_path = Path(f"./Database_DATA/{virus_name}/{pdb_id}.pdb")

            # When PDB, ligand, or chain changes, process the collected functional group info
            if (pdb_id != current_pdb or ligand != current_ligand or chain != current_chain):
                if current_pdb is not None and pdb_file.exists():
                    output_dir = setup_output_directory(base_output_dir, virus_name, current_pdb, current_ligand, current_chain)
                    process_pdb_ligand_chain(current_pdb, current_ligand, current_chain, pdb_file, functional_group_info, output_dir)
                
                # Reset for the new PDB, ligand, and chain
                current_pdb = pdb_id
                current_ligand = ligand
                current_chain = chain
                pdb_file = pdb_path
                functional_group_info = [(functional_group, atom_id, exact_atom, atom_type)]
            else:
                functional_group_info.append((functional_group, atom_id, exact_atom, atom_type))

        # Process the last entry
        if current_pdb is not None and pdb_file.exists():
            output_dir = setup_output_directory(base_output_dir, virus_name, current_pdb, current_ligand, current_chain)
            process_pdb_ligand_chain(current_pdb, current_ligand, current_chain, pdb_file, functional_group_info, output_dir)

# Function to create the PyMOL script
def write_pymol_script(pdb_file, output_script, ligand, chain, functional_groups):
    with open(output_script, 'w') as script:
        # Load the PDB file
        script.write(f"load {pdb_file}\n")
        
        # Create a selection for the ligand
        script.write(f"select {ligand}_Ligand, resn {ligand} and chain {chain}\n")
        script.write(f"create Ligand_Object, {ligand}_Ligand\n")  # Create a new object for the ligand
        script.write(f"show sticks, Ligand_Object\n")  # Show the ligand as sticks
        
        # Create a selection for each functional group, grouping all atoms together
        for fg_name, atoms in functional_groups.items():
            atom_selection = " or ".join([f"id {atom_id}" for atom_id, exact_atom, atom_type in atoms])
            script.write(f"select {fg_name}, chain {chain} and ({atom_selection})\n")
            script.write(f"create {fg_name}_Object, {fg_name}\n")
        
        # Save the PyMOL session
        session_file = Path(output_script).with_suffix('.pse')  # Full path for session file
        script.write(f"save {session_file}\n")
        print(f"PyMOL script created: {output_script}")

# Main entry point
def main():
    base_output_dir = "./Processed_PDBs_FG/"
    virus_names = [
        "Human immunodeficiency virus 1",
        # "Severe acute respiratory syndrome coronavirus 2",
        # "Human papillomavirus",
        # "Human papillomavirus 1",
        # "Human papillomavirus 11",
        # "Human papillomavirus 16",
        # "Human papillomavirus 18",
        # "Human papillomavirus 26",
        # "Human papillomavirus 31",
        # "Human papillomavirus 33",
        # "Human papillomavirus 35",
        # "Human papillomavirus 4",
        # "Human papillomavirus 45",
        # "Human papillomavirus 49",
        # "Human papillomavirus 51",
        # "Human papillomavirus 52",
        # "Human papillomavirus 53",
        # "Human papillomavirus 58",
        # "Human papillomavirus 59",
        # "Human papillomavirus 6",
        # "Human papillomavirus 66",
        # "Human papillomavirus type 16",
        # "Human papillomavirus type 6",
        # "Human papillomavirus type 6a"
     
        # Add more virus names as needed
    ]
    
    for virus_name in virus_names:
        print(f"Processing virus: {virus_name}")
        process_functional_groups_file(virus_name, base_output_dir)

if __name__ == "__main__":
    main()
