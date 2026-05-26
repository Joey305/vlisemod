import csv
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, TimeoutError
import csv
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem, rdFMCS
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import os
import signal
import time
from concurrent.futures import ProcessPoolExecutor, TimeoutError



def find_solvent_exposed_atoms(pdb_file, ligand_resname, ligand_chain, probe_radius=1.40, threshold=0.1):
    # Parse the structure from the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    # Remove water molecules (HOH)
    for model in structure:
        for chain in model:
            residues_to_remove = [residue for residue in chain if residue.resname == "HOH"]
            for residue in residues_to_remove:
                chain.detach_child(residue.id)
    
    # Initialize the Shrake-Rupley algorithm
    sr = ShrakeRupley(probe_radius=probe_radius)
    
    # Compute SASA for the entire structure
    sr.compute(structure, level="A")
    
    # Extract the ligand atoms based on the residue name and chain ID
    ligand_atoms = []
    for model in structure:
        for chain in model:
            if chain.id == ligand_chain:
                for residue in chain:
                    if residue.resname == ligand_resname:
                        ligand_atoms.extend(residue.get_atoms())
    
    # Find solvent-exposed atoms by comparing their SASA values with the threshold
    solvent_exposed_atoms = []
    for atom in ligand_atoms:
        if hasattr(atom, 'sasa') and atom.sasa > threshold:
            solvent_exposed_atoms.append({
                'atom_name': atom.name,
                'atom_id': atom.get_serial_number(),
                'sasa': atom.sasa
            })

    return solvent_exposed_atoms


def save_sasa_data_to_csv(solvent_exposed_atoms, virus_name, pdb_id, ligand_resname, ligand_chain, output_file):
    """
    Saves the solvent-exposed atoms' SASA data into a CSV file, along with Virus Name, PDB ID, Ligand, and Chain.
    """
    file_exists = os.path.isfile(output_file)
    
    with open(output_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        
        # Write header only if the file does not exist
        if not file_exists:
            writer.writerow(["Virus Name", "PDB ID", "Ligand", "Chain", "Atom Name", "Atom ID", "SASA Value"])

        # Write the solvent-exposed atom data
        for atom_info in solvent_exposed_atoms:
            writer.writerow([virus_name, pdb_id, ligand_resname, ligand_chain, atom_info['atom_name'], atom_info['atom_id'], f"{atom_info['sasa']:.2f}"])


def get_batch_data():
    """
    Retrieves virus_name, pdb_id, ligand, and chain from the Ligand_Arp_Diagram table.
    """
    conn = sqlite3.connect("viral_data.db")
    cursor = conn.cursor()
    query = "SELECT virus_name, pdb_id, ligand, chain FROM Ligand_Arp_Diagram"
    cursor.execute(query)
    batch_data = cursor.fetchall()  # Retrieve all data
    conn.close()
    return batch_data


def process_single_entry(entry, output_file):
    """
    Processes a single PDB entry, extracts solvent-exposed atoms, and writes the result to CSV.
    """
    virus_name, pdb_id, ligand_resname, ligand_chain = entry

    try:
        # Skip if this entry already exists in the CSV file
        if os.path.isfile(output_file):
            existing_df = pd.read_csv(output_file, usecols=["Virus Name", "PDB ID", "Ligand", "Chain"])
            if ((existing_df["Virus Name"] == virus_name) & 
                (existing_df["PDB ID"] == pdb_id) &
                (existing_df["Ligand"] == ligand_resname) & 
                (existing_df["Chain"] == ligand_chain)).any():
                print(f"Skipping {virus_name}, {pdb_id}, {ligand_resname}, {ligand_chain} as it already exists.")
                return

        # Find the PDB file location
        pdb_dir = os.path.join("Database_DATA", virus_name)
        pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")
        
        # Check if the PDB file exists
        if not os.path.isfile(pdb_file):
            print(f"PDB file {pdb_file} not found. Skipping.")
            return

        # Find solvent-exposed atoms using Shrake-Rupley algorithm
        solvent_exposed_atoms = find_solvent_exposed_atoms(pdb_file, ligand_resname, ligand_chain)

        # Save results to CSV
        if solvent_exposed_atoms:
            save_sasa_data_to_csv(solvent_exposed_atoms, virus_name, pdb_id, ligand_resname, ligand_chain, output_file)
            print(f"Processed {pdb_id} with ligand {ligand_resname}.")
        else:
            print(f"No solvent-exposed atoms found for {pdb_id} with ligand {ligand_resname}.")

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")


def batch_process_parallel():
    """
    Processes all PDB files in parallel using multiple processes and handles timeouts.
    """
    batch_data = get_batch_data()
    output_file = "Rupley_SASA_data_batch.csv"
    
    # Set the number of processes (e.g., 4 processes)
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(process_single_entry, entry, output_file): entry for entry in batch_data}

        for future in futures:
            try:
                # Set a timeout of 180 seconds (3 minutes) for each task
                future.result(timeout=180)
            except TimeoutError:
                virus_name, pdb_id, ligand_resname, ligand_chain = futures[future]
                print(f"Processing for {virus_name}, {pdb_id}, {ligand_resname}, {ligand_chain} took too long. Skipping this entry.")
            except Exception as e:
                virus_name, pdb_id, ligand_resname, ligand_chain = futures[future]
                print(f"Error processing {pdb_id}: {e}")


if __name__ == "__main__":
    batch_process_parallel()
