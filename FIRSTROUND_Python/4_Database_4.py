import sqlite3
import csv
import os
from pathlib import Path

def initialize_database(db_name="viral_data.db"):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # Create tables with the virus_name column
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ligand_atoms (
            virus_name TEXT,
            pdb_file TEXT,
            ligand TEXT,
            chain TEXT,
            sequence_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ligand_water_distances (
            virus_name TEXT,
            pdb_file TEXT,
            ligand TEXT,
            chain TEXT,
            sequence_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL,
            water_chain TEXT,
            water_sequence_id INTEGER,
            water_x REAL,
            water_y REAL,
            water_z REAL,
            distance REAL
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS solvent_exposed_atoms (
            virus_name TEXT,
            pdb_file TEXT,
            ligand TEXT,
            chain TEXT,
            sequence_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS distal_atoms (
            virus_name TEXT,
            pdb_file TEXT,
            ligand TEXT,
            chain TEXT,
            sequence_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS receptor_binding_pocket (
            virus_name TEXT,
            pdb_file TEXT,
            residue TEXT,
            chain TEXT,
            atom_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL
        );
    ''')

    conn.commit()
    return conn


def adjust_csv_columns(row, table_name):
    if table_name == "ligand_atoms" or table_name == "solvent_exposed_atoms":
        return {
            'virus_name': row.get('virus_name', None),
            'pdb_file': row.get('PDB File', None),
            'ligand': row.get('Ligand', row.get('Residue', None)),
            'chain': row.get('Chain', None),
            'sequence_id': row.get('Sequence ID', row.get('Atom ID', None)),
            'exact_atom': row.get('Exact Atom', None),
            'atom_type': row.get('Atom Type', None),
            'x': row.get('X', None),
            'y': row.get('Y', None),
            'z': row.get('Z', None)
        }
    elif table_name == "ligand_water_distances":
        return {
            'virus_name': row.get('virus_name', None),
            'pdb_file': row.get('PDB File', None),
            'ligand': row.get('Ligand', row.get('Residue', None)),
            'chain': row.get('Chain', None),
            'sequence_id': row.get('Sequence ID', row.get('Atom ID', None)),
            'exact_atom': row.get('Exact Atom', None),
            'atom_type': row.get('Atom Type', None),
            'x': row.get('X', None),
            'y': row.get('Y', None),
            'z': row.get('Z', None),
            'water_chain': row.get('Water Chain', None),
            'water_sequence_id': row.get('Water Sequence ID', None),
            'water_x': row.get('Water X', None),
            'water_y': row.get('Water Y', None),
            'water_z': row.get('Water Z', None),
            'distance': row.get('Distance', None)
        }
    elif table_name == "receptor_binding_pocket":
        return {
            'virus_name': row.get('virus_name', None),
            'pdb_file': row.get('PDB File', None),
            'residue': row.get('Residue', None),
            'chain': row.get('Chain', None),
            'atom_id': row.get('Atom ID', None),
            'exact_atom': row.get('Exact Atom', None),
            'atom_type': row.get('Atom Type', None),
            'x': row.get('X', None),
            'y': row.get('Y', None),
            'z': row.get('Z', None)
        }


def insert_csv_data(conn, csv_file, table_name, virus_name):
    cursor = conn.cursor()
    try:
        with open(csv_file, 'r') as file:
            # Check if the file has headers
            first_line = file.readline().strip().split(',')
            expected_headers = ["PDB File", "Ligand", "Chain", "Sequence ID", "Exact Atom", "Atom Type", "X", "Y", "Z"]

            if first_line != expected_headers:
                print(f"Warning: Missing or incorrect headers in file {csv_file}. Skipping...")
                return

            # Move file pointer back to start after checking headers
            file.seek(0)

            reader = csv.DictReader(file)
            if reader.fieldnames is None:
                print(f"Warning: No data found in file {csv_file}. Skipping...")
                return

            rows = []
            for row in reader:
                row['virus_name'] = virus_name
                adjusted_row = adjust_csv_columns(row, table_name)
                rows.append(tuple(adjusted_row.values()))

            if not rows:
                print(f"Warning: No rows found in file {csv_file}. Skipping...")
                return

            placeholders = ', '.join(['?'] * len(rows[0]))
            query = f"INSERT INTO {table_name} VALUES ({placeholders})"
            cursor.executemany(query, rows)
            print(f"Data successfully inserted from {csv_file} into {table_name}")
    except Exception as e:
        print(f"Error processing file {csv_file}: {e}")
    conn.commit()


def load_data_into_database(base_directory, db_name="viral_data.db"):
    conn = initialize_database(db_name)
    base_path = Path(base_directory)

    for virus_dir in base_path.iterdir():
        if virus_dir.is_dir():
            virus_name = virus_dir.name  # Capture the virus name
            print(f"Processing virus: {virus_name}")
            for ligand_dir in virus_dir.iterdir():
                if ligand_dir.is_dir():
                    print(f"Processing ligand: {ligand_dir.name}")
                    for csv_file in ligand_dir.glob("*.csv"):
                        print(f"Processing file: {csv_file}")
                        if "ligand.csv" in csv_file.name:
                            insert_csv_data(conn, csv_file, "ligand_atoms", virus_name)
                        elif "ligand_water_distances.csv" in csv_file.name:
                            insert_csv_data(conn, csv_file, "ligand_water_distances", virus_name)
                        elif "solvent_exposed_atoms.csv" in csv_file.name:
                            insert_csv_data(conn, csv_file, "solvent_exposed_atoms", virus_name)
                        elif "distalAtoms.csv" in csv_file.name:
                            insert_csv_data(conn, csv_file, "distal_atoms", virus_name)
                        elif "receptor_binding_pocket.csv" in csv_file.name:
                            insert_csv_data(conn, csv_file, "receptor_binding_pocket", virus_name)

    print(f"All data loaded into {db_name}.")
    conn.close()


# Example usage
base_directory = r"C:\Users\joeys\OneDrive\Documents\Projects\VIRAL_DATABASE\Processed_PDBs"
load_data_into_database(base_directory)
