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




class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException




def read_ligand_expo():
    file_name = "Components-smiles-stereo-oe.smi"
    try:
        df = pd.read_csv(file_name, sep="\t", header=None, names=["SMILES", "ID", "Name"])
        df.set_index("ID", inplace=True)
        return df.to_dict()["SMILES"]
    except FileNotFoundError:
        print(f"File {file_name} not found. Please ensure the file is in the correct directory.")
        return {}


def generate_conformer(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Failed to generate RDKit molecule from SMILES: {smiles_string}")
        return None

    mol = Chem.AddHs(mol)
    
    # Check the result of the embedding process
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:  # If result is not 0, embedding failed
        print(f"Embedding failed for SMILES: {smiles_string}")
        return None

    try:
        AllChem.UFFOptimizeMolecule(mol)
    except ValueError as e:
        print(f"Error during UFF optimization for SMILES: {smiles_string}, error: {e}")
        return None

    mol = Chem.RemoveHs(mol)  # Remove hydrogens after conformer generation
    return mol


def match_atoms_using_mcs(pdb_mol, smiles_mol, atom_mapping):
    """
    Match atoms from the PDB molecule and SMILES conformer using Maximum Common Substructure (MCS).
    """
    # Use RDKit's MCS to find the best match between the PDB molecule and SMILES conformer
    res = rdFMCS.FindMCS([pdb_mol, smiles_mol], completeRingsOnly=True, bondCompare=rdFMCS.BondCompare.CompareOrder)
    patt = Chem.MolFromSmarts(res.smartsString)
    
    if patt is None:
        print("MCS pattern could not be generated.")
        return []
    
    pdb_matches = pdb_mol.GetSubstructMatch(patt)
    smiles_matches = smiles_mol.GetSubstructMatch(patt)

    if len(pdb_matches) == 0 or len(smiles_matches) == 0:
        print("No matches found between PDB and SMILES conformer.")
        return []

    # Use the atom_mapping to retain the correct PDB atom names and IDs
    matches = [(atom_mapping[pdb_idx][0], atom_mapping[pdb_idx][1], pdb_idx, smiles_idx) for pdb_idx, smiles_idx in zip(pdb_matches, smiles_matches)]
    
    print(f"Found {len(matches)} atom matches using MCS.")
    return matches, pdb_matches, smiles_matches


def parse_pdb_file(pdb_file, ligand_resname):
    """
    Parse PDB file and extract ligand atoms with atom name and atom ID.
    """
    pdb_mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if pdb_mol is None:
        print(f"Failed to parse PDB file: {pdb_file}")
        return None, None

    # Map PDB atom index to (PDB atom name, PDB atom ID)
    atom_mapping = {}
    for atom in pdb_mol.GetAtoms():
        res_name = atom.GetPDBResidueInfo().GetResidueName().strip()
        if res_name == ligand_resname:
            pdb_atom_name = atom.GetPDBResidueInfo().GetName().strip()
            pdb_atom_id = atom.GetPDBResidueInfo().GetSerialNumber()
            atom_mapping[atom.GetIdx()] = (pdb_atom_name, pdb_atom_id)

    return pdb_mol, atom_mapping


def save_to_csv(matches, output_file, virus_name, pdb_id, ligand, chain):
    """
    Appends the matched atom data into a CSV file with additional columns for virus_name, pdb_id, ligand, and chain.
    If the file does not exist, it writes the header first.
    """
    # Check if the file already exists to avoid rewriting the header
    file_exists = os.path.isfile(output_file)

    with open(output_file, mode='a', newline='') as file:  # Open in append mode
        writer = csv.writer(file)

        # Write the header if the file doesn't exist
        if not file_exists:
            writer.writerow(["Virus Name", "PDB ID", "Ligand", "Chain", "PDB Atom Name", "PDB Atom ID", "PDB Atom Index", "SMILES Atom Index"])

        # Write the data for each match, including the virus_name, pdb_id, ligand, and chain
        for match in matches:
            writer.writerow([virus_name, pdb_id, ligand, chain, *match])

def save_unmatched_atoms_to_csv(unmatched_pdb, unmatched_smiles, pdb_mapping, smiles_mol, output_file):
    """
    Saves the unmatched PDB atoms and SMILES atoms into a CSV file.
    """
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Unmatched PDB Atom Name", "Unmatched PDB Atom ID", "Unmatched SMILES Atom Index"])

        # Write unmatched PDB atoms
        for pdb_idx in unmatched_pdb:
            pdb_atom_name, pdb_atom_id = pdb_mapping[pdb_idx]
            writer.writerow([pdb_atom_name, pdb_atom_id, "N/A"])

        # Write unmatched SMILES atoms
        for smiles_idx in unmatched_smiles:
            writer.writerow(["N/A", "N/A", smiles_idx])


def find_bonded_carbons(pdb_mol, atom_mapping, unmatched_pdb_indices):
    """
    For each unmatched PDB oxygen, find the carbon atom(s) it is bonded to.
    """
    bonded_carbons = {}
    for pdb_idx in unmatched_pdb_indices:
        pdb_atom = pdb_mol.GetAtomWithIdx(pdb_idx)
        if pdb_atom.GetSymbol() == "O":  # Only handle oxygen atoms
            for neighbor in pdb_atom.GetNeighbors():
                if neighbor.GetSymbol() == "C":  # Find carbon atoms
                    carbon_idx = neighbor.GetIdx()
                    bonded_carbons[pdb_idx] = carbon_idx
    return bonded_carbons


def find_smiles_carbon_oxygen(smiles_mol, carbon_idx):
    """
    For a given carbon atom in the SMILES molecule, find the double-bonded oxygen atom.
    """
    carbon_atom = smiles_mol.GetAtomWithIdx(carbon_idx)
    for bond in carbon_atom.GetBonds():
        neighbor = bond.GetOtherAtom(carbon_atom)
        if neighbor.GetSymbol() == "O" and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            return neighbor.GetIdx()  # Return the SMILES oxygen index
    return None


def match_unmatched_oxygens(pdb_mol, smiles_mol, atom_mapping, unmatched_pdb_indices, unmatched_smiles_indices, matched_atoms):
    """
    Find the corresponding oxygen atom in SMILES for the unmatched PDB oxygen atoms by checking their bonded carbons.
    """
    # Step 1: Find the bonded carbon atoms for the unmatched PDB oxygens
    bonded_carbons = find_bonded_carbons(pdb_mol, atom_mapping, unmatched_pdb_indices)

    new_matches = []
    
    # Step 2: For each bonded carbon in PDB, find the corresponding SMILES carbon
    for pdb_oxygen_idx, pdb_carbon_idx in bonded_carbons.items():
        # Check if this PDB carbon atom is already matched to a SMILES carbon atom
        matched_smiles_carbon_idx = None
        for match in matched_atoms:
            pdb_atom_name, pdb_atom_id, pdb_idx, smiles_idx = match
            if pdb_idx == pdb_carbon_idx:
                matched_smiles_carbon_idx = smiles_idx
                break

        if matched_smiles_carbon_idx is not None:
            # Step 3: Find the double-bonded oxygen in SMILES
            matched_smiles_oxygen_idx = find_smiles_carbon_oxygen(smiles_mol, matched_smiles_carbon_idx)
            if matched_smiles_oxygen_idx is not None and matched_smiles_oxygen_idx in unmatched_smiles_indices:
                # Add the new match to the list
                pdb_atom_name = pdb_mol.GetAtomWithIdx(pdb_oxygen_idx).GetPDBResidueInfo().GetName().strip()
                pdb_atom_id = pdb_mol.GetAtomWithIdx(pdb_oxygen_idx).GetPDBResidueInfo().GetSerialNumber()
                new_matches.append((pdb_atom_name, pdb_atom_id, pdb_oxygen_idx, matched_smiles_oxygen_idx))

    return new_matches


def get_batch_data():
    """
    Retrieves virus_name, pdb_id, ligand, and chain from the Ligand_Arp_Diagram table
    """
    conn = sqlite3.connect("viral_data.db")
    cursor = conn.cursor()
    query = "SELECT virus_name, pdb_id, ligand, chain FROM Ligand_Arp_Diagram"
    cursor.execute(query)
    batch_data = cursor.fetchall()  # Retrieve all data
    conn.close()
    return batch_data

def get_smiles_from_db(virus_name, pdb_id, ligand):
    """
    Retrieve SMILES from the Functional_GROUPED table based on virus_name, pdb_id, and ligand
    """
    conn = sqlite3.connect("viral_data.db")
    cursor = conn.cursor()
    query = """SELECT smiles FROM Functional_GROUPED WHERE virus_name=? AND pdb_file=? AND ligand=?"""
    cursor.execute(query, (virus_name, pdb_id, ligand))
    result = cursor.fetchone()
    conn.close()
    return result[0] if result else None



import os
import pandas as pd

def save_to_csv_if_new(matches, output_file, virus_name, pdb_id, ligand, chain):
    """
    Appends matched atom data into a CSV file only if a row with the same virus_name, pdb_id, ligand, and chain
    does not already exist in the CSV.
    """
    # Check if the file already exists
    file_exists = os.path.isfile(output_file)

    # If the file exists, load only the relevant columns: virus_name, pdb_id, ligand, and chain
    if file_exists:
        existing_df = pd.read_csv(output_file, usecols=["Virus Name", "PDB ID", "Ligand", "Chain"])
        # Check if the combination of virus_name, pdb_id, ligand, and chain already exists
        if ((existing_df["Virus Name"] == virus_name) &
            (existing_df["PDB ID"] == pdb_id) &
            (existing_df["Ligand"] == ligand) &
            (existing_df["Chain"] == chain)).any():
            print(f"Skipping {virus_name}, {pdb_id}, {ligand}, {chain} as it already exists in the file.")
            return  # Skip the current entry if it already exists
    else:
        # If file doesn't exist, initialize an empty DataFrame for checking
        existing_df = pd.DataFrame(columns=["Virus Name", "PDB ID", "Ligand", "Chain"])

    # Create a DataFrame for new matches
    new_data = []
    for match in matches:
        new_data.append([virus_name, pdb_id, ligand, chain, *match])

    new_df = pd.DataFrame(new_data, columns=["Virus Name", "PDB ID", "Ligand", "Chain", "PDB Atom Name", "PDB Atom ID", "PDB Atom Index", "SMILES Atom Index"])

    # Append new data to the CSV
    new_df.to_csv(output_file, mode='a', header=not file_exists, index=False)
    print(f"Processed and saved {len(new_data)} new rows to {output_file}.")

def batch_process():
    """
    This function processes each PDB file and ligand in batch, ensuring that no duplicates are saved.
    The duplicate check happens before any computational steps like generating conformers.
    If processing takes longer than 3 minutes, the file is skipped.
    """
    batch_data = get_batch_data()  # Fetch batch data from database
    output_file = "PDB_SMILE_Code_Batch.csv"  # The output CSV file

    # Set up signal for timeout (180 seconds)
    signal.signal(signal.SIGALRM, timeout_handler)

    for virus_name, pdb_id, ligand, chain in batch_data:
        try:
            # Start the timer for 180 seconds (3 minutes)
            signal.alarm(180)

            # Step 1: Check if this combination of virus_name, pdb_id, ligand, and chain already exists in the output CSV
            file_exists = os.path.isfile(output_file)
            if file_exists:
                # Load only relevant columns for checking duplicates
                existing_df = pd.read_csv(output_file, usecols=["Virus Name", "PDB ID", "Ligand", "Chain"])
                if ((existing_df["Virus Name"] == virus_name) &
                    (existing_df["PDB ID"] == pdb_id) &
                    (existing_df["Ligand"] == ligand) &
                    (existing_df["Chain"] == chain)).any():
                    print(f"Skipping {virus_name}, {pdb_id}, {ligand}, {chain} as it already exists in the file.")
                    continue  # Skip this entry if it already exists

            # Step 2: Proceed with processing only if the row does not already exist
            pdb_dir = os.path.join("Database_DATA", virus_name)  # Directory containing PDB files
            pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")  # Full path to the PDB file

            # Get SMILES from the database
            smiles_string = get_smiles_from_db(virus_name, pdb_id, ligand)
            if smiles_string is None:
                print(f"No SMILES found for {ligand} in {pdb_id}, skipping...")
                continue

            # Generate 3D conformer from SMILES
            smiles_mol = generate_conformer(smiles_string)
            if smiles_mol is None:
                print(f"Failed to generate conformer for {ligand} in {pdb_id}, skipping...")
                continue

            # Parse PDB file and match atoms using MCS
            pdb_mol, atom_mapping = parse_pdb_file(pdb_file, ligand)
            if pdb_mol is None:
                print(f"Failed to parse PDB: {pdb_file}, skipping...")
                continue

            # Initial MCS atom matching
            matched_atoms, pdb_matches, smiles_matches = match_atoms_using_mcs(pdb_mol, smiles_mol, atom_mapping)

            # Unmatched PDB and SMILES atoms
            unmatched_pdb_indices = set(atom_mapping.keys()) - set(pdb_matches)
            unmatched_smiles_indices = set(range(smiles_mol.GetNumAtoms())) - set(smiles_matches)

            # Match unmatched oxygen atoms
            new_oxygen_matches = match_unmatched_oxygens(pdb_mol, smiles_mol, atom_mapping, unmatched_pdb_indices, unmatched_smiles_indices, matched_atoms)

            # Append the new oxygen matches to the existing matches
            matched_atoms.extend(new_oxygen_matches)

            # Save results to CSV using the new function that skips duplicates
            save_to_csv_if_new(matched_atoms, output_file, virus_name, pdb_id, ligand, chain)
            print(f"Processed {pdb_id} with ligand {ligand} successfully.")

            # Disable the alarm after successful processing
            signal.alarm(0)

        except TimeoutException:
            print(f"Processing for {virus_name}, {pdb_id}, {ligand}, {chain} took too long. Skipping this entry.")
            signal.alarm(0)  # Disable the alarm
            continue

        except Exception as e:
            print(f"Error processing {pdb_id}: {e}")
            continue

def process_single_entry(entry):
    """
    Process a single batch entry (virus_name, pdb_id, ligand, chain).
    This function will be called in parallel for different entries.
    """
    virus_name, pdb_id, ligand, chain = entry
    output_file = "PDB_SMILE_Code_Batch.csv"

    try:
        file_exists = os.path.isfile(output_file)
        if file_exists:
            existing_df = pd.read_csv(output_file, usecols=["Virus Name", "PDB ID", "Ligand", "Chain"])
            if ((existing_df["Virus Name"] == virus_name) &
                (existing_df["PDB ID"] == pdb_id) &
                (existing_df["Ligand"] == ligand) &
                (existing_df["Chain"] == chain)).any():
                print(f"Skipping {virus_name}, {pdb_id}, {ligand}, {chain} as it already exists in the file.")
                return

        pdb_dir = os.path.join("Database_DATA", virus_name)
        pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")

        smiles_string = get_smiles_from_db(virus_name, pdb_id, ligand)
        if smiles_string is None:
            print(f"No SMILES found for {ligand} in {pdb_id}, skipping...")
            return

        smiles_mol = generate_conformer(smiles_string)
        if smiles_mol is None:
            print(f"Failed to generate conformer for {ligand} in {pdb_id}, skipping...")
            return

        pdb_mol, atom_mapping = parse_pdb_file(pdb_file, ligand)
        if pdb_mol is None:
            print(f"Failed to parse PDB: {pdb_file}, skipping...")
            return

        matched_atoms, pdb_matches, smiles_matches = match_atoms_using_mcs(pdb_mol, smiles_mol, atom_mapping)
        unmatched_pdb_indices = set(atom_mapping.keys()) - set(pdb_matches)
        unmatched_smiles_indices = set(range(smiles_mol.GetNumAtoms())) - set(smiles_matches)
        new_oxygen_matches = match_unmatched_oxygens(pdb_mol, smiles_mol, atom_mapping, unmatched_pdb_indices, unmatched_smiles_indices, matched_atoms)

        matched_atoms.extend(new_oxygen_matches)
        save_to_csv_if_new(matched_atoms, output_file, virus_name, pdb_id, ligand, chain)
        print(f"Processed {pdb_id} with ligand {ligand} successfully.")

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")


def batch_process_parallel():
    """
    Processes all PDB files in parallel using multiple processes and handles timeouts.
    """
    batch_data = get_batch_data()
    
    # Set the number of processes (e.g., 4 processes)
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(process_single_entry, entry): entry for entry in batch_data}

        for future in futures:
            try:
                # Set a timeout of 180 seconds (3 minutes) for each task
                future.result(timeout=180)
            except TimeoutError:
                virus_name, pdb_id, ligand, chain = futures[future]
                print(f"Processing for {virus_name}, {pdb_id}, {ligand}, {chain} took too long. Skipping this entry.")
            except Exception as e:
                virus_name, pdb_id, ligand, chain = futures[future]
                print(f"Error processing {pdb_id}: {e}")


if __name__ == "__main__":
    batch_process_parallel()