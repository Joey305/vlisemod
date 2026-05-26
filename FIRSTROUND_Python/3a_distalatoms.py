import sqlite3
import csv
from pathlib import Path

def initialize_database(db_name="viral_data.db"):
    # Ensure the distal_atoms table exists
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

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

    conn.commit()
    return conn

def insert_distal_atoms_data(conn, csv_file, virus_name):
    cursor = conn.cursor()
    try:
        with open(csv_file, 'r') as file:
            reader = csv.DictReader(file)

            # Ensure the correct headers are present
            expected_headers = ['PDB File', 'Ligand', 'Chain', 'Sequence ID', 'Exact Atom', 'Atom Type', 'X', 'Y', 'Z']
            if reader.fieldnames is None or not all(header in reader.fieldnames for header in expected_headers):
                print(f"Warning: Missing or incorrect headers in file {csv_file}. Skipping...")
                return

            rows = []
            for row in reader:
                # Ensure that 10 values are being provided (including virus_name)
                rows.append((
                    virus_name,             # Dynamically add virus_name
                    row['PDB File'],        # From CSV
                    row['Ligand'],          # From CSV
                    row['Chain'],           # From CSV
                    row['Sequence ID'],     # From CSV
                    row['Exact Atom'],      # From CSV
                    row['Atom Type'],       # From CSV
                    row['X'],               # From CSV
                    row['Y'],               # From CSV
                    row['Z']                # From CSV
                ))

            if not rows:
                print(f"Warning: No rows found in file {csv_file}. Skipping...")
                return

            placeholders = ', '.join(['?'] * 10)  # Now 10 placeholders for 10 columns
            query = f"INSERT INTO distal_atoms VALUES ({placeholders})"
            cursor.executemany(query, rows)
            print(f"Data successfully inserted from {csv_file} into distal_atoms")
    except Exception as e:
        print(f"Error processing file {csv_file}: {e}")
    conn.commit()

# Traverse the directory structure and load distal_atoms CSV data
def load_distal_atoms_data(base_directory, db_name="viral_data.db"):
    conn = initialize_database(db_name)
    base_path = Path(base_directory)

    for virus_dir in base_path.iterdir():
        if virus_dir.is_dir():
            virus_name = virus_dir.name  # Capture the virus name
            print(f"Processing virus: {virus_name}")
            for ligand_dir in virus_dir.iterdir():
                if ligand_dir.is_dir():
                    print(f"Processing ligand: {ligand_dir.name}")
                    for csv_file in ligand_dir.glob("*distalAtoms.csv"):
                        print(f"Processing file: {csv_file}")
                        insert_distal_atoms_data(conn, csv_file, virus_name)

    print(f"All distal_atoms data loaded into {db_name}.")
    conn.close()

# Example usage
base_directory = r"C:\Users\joeys\OneDrive\Documents\Projects\VIRAL_DATABASE\Processed_PDBs"
load_distal_atoms_data(base_directory)
