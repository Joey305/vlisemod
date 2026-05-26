import os
import csv
import logging
from pathlib import Path
from math import sqrt

# Configure logging
logging.basicConfig(
    filename='Add_Carbonyl_log.log',
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.DEBUG
)

# Function to calculate distance between two points in 3D space
def calculate_distance(coords1, coords2):
    return sqrt((coords1[0] - coords2[0]) ** 2 + (coords1[1] - coords2[1]) ** 2 + (coords1[2] - coords2[2]) ** 2)

# Function to find the closest missing oxygen atom using both adjacency and 3D distance
def find_closest_oxygen(pdb_id, ligand_id, chain, group_name, atom_ids, ligand_atoms, appended_atoms, functional_groups):
    closest_oxygen_id = None
    min_avg_distance = float('inf')

    # First, try to find oxygen atoms within the +/- 2 range
    for atom_id in atom_ids:
        for offset in [-2, -1, 1, 2]:  # Adjust as necessary to find adjacent atoms
            adjacent_atom_id = atom_id + offset
            
            # Skip if the adjacent atom is already part of the functional group
            if adjacent_atom_id in atom_ids:
                continue
            
            adjacent_atom = ligand_atoms.get((pdb_id, ligand_id, chain, adjacent_atom_id))
            
            if adjacent_atom and adjacent_atom_id not in appended_atoms:
                if adjacent_atom[1] == 'O':  # Check if the atom is an Oxygen
                    oxygen_coords = [float(adjacent_atom[i]) for i in range(2, 5)]
                    distances = []

                    for group_atom_id in atom_ids:
                        if (pdb_id, ligand_id, chain, group_atom_id) not in ligand_atoms:
                            logging.warning(f"Atom ID {group_atom_id} not found in ligand_atoms for PDB ID: {pdb_id}, Ligand ID: {ligand_id}, Chain: {chain}")
                            continue
                        
                        group_atom_coords = ligand_atoms.get((pdb_id, ligand_id, chain, group_atom_id))
                        if group_atom_coords is None:
                            logging.warning(f"Group atom ID {group_atom_id} not found in ligand_atoms for PDB ID: {pdb_id}, Ligand ID: {ligand_id}, Chain: {chain}")
                            continue

                        group_atom_coords = [float(group_atom_coords[i]) for i in range(2, 5)]
                        distances.append(calculate_distance(oxygen_coords, group_atom_coords))
                    
                    if distances:  # Ensure there are distances to average
                        avg_distance = sum(distances) / len(distances)

                        if avg_distance < min_avg_distance:
                            min_avg_distance = avg_distance
                            closest_oxygen_id = adjacent_atom_id

    # If no oxygen atom was found within +/- 2 range, use 3D distance for all oxygen atoms
    if closest_oxygen_id is None:
        for atom_id in atom_ids:
            group_atom_coords = ligand_atoms.get((pdb_id, ligand_id, chain, atom_id))
            if group_atom_coords is None:
                logging.warning(f"Group atom ID {atom_id} not found in ligand_atoms for PDB ID: {pdb_id}, Ligand ID: {ligand_id}, Chain: {chain}")
                continue

            group_atom_coords = [float(group_atom_coords[i]) for i in range(2, 5)]
            for key, atom_data in ligand_atoms.items():
                if key[0] == pdb_id and key[1] == ligand_id and key[2] == chain and atom_data[1] == 'O':  # Check if it’s an oxygen atom
                    oxygen_coords = [float(atom_data[i]) for i in range(2, 5)]
                    distance = calculate_distance(group_atom_coords, oxygen_coords)

                    if distance < min_avg_distance:
                        min_avg_distance = distance
                        closest_oxygen_id = key[3]  # Get the sequence_id of the closest oxygen atom

    return closest_oxygen_id

# Function to append "Exact Atom" and missing oxygen to the output
def append_exact_atom_and_missing_oxygen(input_file, output_file, ligand_atoms_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(ligand_atoms_file, 'r') as ligand_file:
        reader = csv.reader(infile, delimiter='\t')
        ligand_reader = csv.reader(ligand_file)

        # Skip the header rows
        next(reader)
        next(ligand_reader)

        # Load the ligand atom data from ligand_atoms.csv
        ligand_atoms = {
            (row[1], row[2], row[3], int(row[4])): (row[5], row[6], row[7], row[8], row[9])  # (pdb_file, ligand, chain, sequence_id): (exact_atom, atom_type, x, y, z)
            for row in ligand_reader
        }

        # Write the header to the output text file
        header = "PDB ID\tLigand ID\tChain\tFunctional Group\tAtom ID\tExact Atom\tAtom Type\n"
        outfile.write(header)

        # Create a temporary storage for functional groups
        functional_groups = {}
        appended_atoms_global = set()
        last_ligand_id = None  # To track ligand changes

        # Define expected atom counts for common functional groups
        expected_atom_counts = {
        'Hydroxyl/Alcohol': 1,  # 1 oxygen
        'Ether': 3,  # 1 oxygen + 2 carbon
        'Aldehyde': 3,  # 1 carbon + 1 oxygen + 1 hydrogen
        'Ketone': 3,  # 1 carbon + 1 oxygen + 1 carbon
        'Carboxylic-Acid': 3,  # 1 carbon + 2 oxygen
        'Ester': 4,  # 1 carbon + 2 oxygen + 1 carbon
        'Terminal-Methyl-Ester': 5,  # 1 carbon + 2 oxygen + 1 carbon + 1 hydrogen
        'Amine': 1,  # 1 nitrogen
        'Amide': 4,  # 1 nitrogen + 1 carbon + 1 oxygen + 1 carbon
        'Benzene': 6,  # 6 carbon
        'Sulfonamide': 4,  # 1 sulfur + 2 oxygen + 1 nitrogen
        'Thiol': 1,  # 1 sulfur
        'Nitrile': 2,  # 1 carbon + 1 nitrogen
        'Nitro': 3,  # 1 nitrogen + 2 oxygen
        'Azo': 2,  # 2 nitrogen
        'Hydrazine': 2,  # 2 nitrogen
        'Isonitrile': 2,  # 1 carbon + 1 nitrogen
        'Imide': 5,  # 2 carbon + 2 oxygen + 1 nitrogen
        'Sulfoxide': 2,  # 1 sulfur + 1 oxygen
        'Sulfone': 3,  # 1 sulfur + 2 oxygen
        'Phosphoric-Acid': 4,  # 1 phosphorus + 3 oxygen
        'Sulfenic-Acid': 2,  # 1 sulfur + 1 oxygen
        'Terminal-Alkyne': 2,  # 2 carbon
        'Aryl-Halide': 7,  # 6 carbon + 1 halide (F, Cl, Br, I)
        'Terminal-Alkene': 2,  # 2 carbon
        'Thioester': 3,  # 1 carbon + 1 sulfur + 1 oxygen
        'Isocyanate': 3,  # 1 nitrogen + 1 carbon + 1 oxygen
        'Aldoxime': 4,  # 1 carbon + 1 oxygen + 1 nitrogen + 1 hydrogen
        'Carbamate': 4,  # 1 nitrogen + 1 carbon + 2 oxygen
        'Isothiocyanate': 3,  # 1 nitrogen + 1 carbon + 1 sulfur
        'Phosphonate-Terminal': 3,  # 1 phosphorus + 2 oxygen
        'Haloalkane': 2,  # 1 carbon + 1 halide (F, Cl, Br, I)
        'Diazo': 2,  # 2 nitrogen
        'Azide': 3,  # 3 nitrogen
        'Phenol': 7,  # 6 carbon + 1 oxygen
        'Thioether': 3,  # 2 carbon + 1 sulfur
        'Carbamic-Acid': 4,  # 1 nitrogen + 1 carbon + 2 oxygen
        'Urea': 4,  # 2 nitrogen + 1 carbon + 1 oxygen
        'Amidine': 4,  # 2 nitrogen + 1 carbon + 1 hydrogen
        'Guanidine': 5,  # 3 nitrogen + 1 carbon + 1 hydrogen
        'Sulfonic-Acid': 4,  # 1 sulfur + 3 oxygen
        'Thiol-Ester': 4,  # 2 carbon + 1 sulfur + 1 oxygen
        'Phosphine': 1,  # 1 phosphorus
        'Imine': 2,  # 1 nitrogen + 1 carbon
        'Hydrazone': 3,  # 2 nitrogen + 1 carbon
        'Sulfonic-Acid-Ester': 5,  # 1 sulfur + 3 oxygen + 1 carbon
        'Thiocyanate': 3,  # 1 sulfur + 1 carbon + 1 nitrogen
        'Carbodiimide': 3,  # 2 nitrogen + 1 carbon
        'Selenide': 2,  # 1 selenium + 1 carbon
        'Acyl-Halide': 3,  # 1 carbon + 1 oxygen + 1 halide (F, Cl, Br, I)
        'Enamine': 3,  # 1 nitrogen + 2 carbon
        'Boronic-Acid': 4,  # 1 boron + 3 oxygen
        'Pyridine': 6,  # 5 carbon + 1 nitrogen
        'Pyrrole': 5,  # 4 carbon + 1 nitrogen
        'Oxime': 3,  # 1 carbon + 1 nitrogen + 1 oxygen
        'Cyclohexane': 6,  # 6 carbon atoms in the ring
        }





        for row in reader:
            if row:  # Make sure the row is not empty
                pdb_id = row[0]
                ligand_id = row[1]
                chain = row[2]
                group_name = row[3]
                atom_id = int(row[4])  # This should now correctly match to PDB Atom ID

                # Store the information by functional group
                group_key = (pdb_id, ligand_id, chain, group_name)
                if group_key not in functional_groups:
                    functional_groups[group_key] = []
                functional_groups[group_key].append(atom_id)

                # Identify the exact atom name and atom type
                exact_atom, atom_type = None, None
                if (pdb_id, ligand_id, chain, atom_id) in ligand_atoms:
                    exact_atom, atom_type = ligand_atoms[(pdb_id, ligand_id, chain, atom_id)][:2]  # Unpack only the first two values
                else:
                    logging.warning(f"No match found for PDB ID: {pdb_id}, Atom ID: {atom_id} in ligand_atoms.csv.")

                # Write the row to the output file
                outfile.write(f"{pdb_id}\t{ligand_id}\t{chain}\t{group_name}\t{atom_id}\t{exact_atom}\t{atom_type}\n")

        # After reading all rows, process the functional groups to find and add missing atoms
        for group_key, atom_ids in functional_groups.items():
            pdb_id, ligand_id, chain, group_name = group_key

            appended_atoms = set()  # Set to keep track of already appended atoms

            base_group_name = group_name.split(' ')[0]
            actual_count = len(atom_ids)

            if base_group_name in expected_atom_counts:
                expected_atoms = expected_atom_counts[base_group_name]
                logging.info(f"Functional group '{group_name}' in PDB ID: {pdb_id}, Ligand ID: {ligand_id}, Chain: {chain} has {actual_count} atoms; Expected: {expected_atoms}")

                if actual_count < expected_atoms:
                    logging.warning(f"Functional group '{group_name}' in PDB ID: {pdb_id}, Ligand ID: {ligand_id}, Chain: {chain} is missing atoms.")
                    missing_atom_count = expected_atoms - actual_count
                    
                    for _ in range(missing_atom_count):  # Try to append the missing oxygen atoms
                        closest_oxygen_id = find_closest_oxygen(pdb_id, ligand_id, chain, group_name, atom_ids, ligand_atoms, appended_atoms, functional_groups)
                        if closest_oxygen_id and closest_oxygen_id not in atom_ids:
                            exact_atom, atom_type = ligand_atoms.get((pdb_id, ligand_id, chain, closest_oxygen_id), (None, None))[:2]  # Unpack only the first two values
                            if (group_key, closest_oxygen_id) not in appended_atoms_global:
                                outfile.write(f"{pdb_id}\t{ligand_id}\t{chain}\t{group_name}\t{closest_oxygen_id}\t{exact_atom}\t{atom_type}\n")
                                logging.info(f"Appended missing atom {closest_oxygen_id} for functional group '{group_name}' in ligand with PDB ID: {pdb_id}, Ligand ID: {ligand_id}, Chain: {chain}")
                                appended_atoms_global.add((group_key, closest_oxygen_id))
                                atom_ids.append(closest_oxygen_id)  # Add to the atom_ids to ensure it's not reappended

# Main function
def main():
    base_input_dir = './Database_DATA'
    ligand_atoms_file = 'output_csvs/ligand_atoms.csv'  # Updated path to your ligand atoms file

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
        # Add more virus names as needed
    ]

    # Process each virus and append missing atoms
    for virus_name in virus_names:
        smiles_file = Path(base_input_dir) / 'Sorted_PDB_Files' / virus_name / 'with_ligands_smiles_functional_groups.txt'
        output_file = smiles_file.with_name(smiles_file.stem + '-appended.txt')
        
        if not smiles_file.exists():
            logging.error(f"SMILES file not found for virus '{virus_name}': {smiles_file}")
            continue

        logging.info(f"Processing SMILES file for virus '{virus_name}': {smiles_file}")
        append_exact_atom_and_missing_oxygen(str(smiles_file), str(output_file), ligand_atoms_file)

if __name__ == "__main__":
    main()

