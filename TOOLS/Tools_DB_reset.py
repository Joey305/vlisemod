import sqlite3
import csv
import os
from pathlib import Path

# Initialize the SQLite database with corrected schema
def initialize_database_with_virus_name(db_name="viral_data.db"):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # Drop existing tables if they exist
    cursor.execute('DROP TABLE IF EXISTS ligand_atoms')
    cursor.execute('DROP TABLE IF EXISTS ligand_water_distances')
    cursor.execute('DROP TABLE IF EXISTS solvent_exposed_atoms')
    cursor.execute('DROP TABLE IF EXISTS distal_atoms')
    cursor.execute('DROP TABLE IF EXISTS receptor_binding_pocket')

    # Recreate the tables with the correct column names and virus_name
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ligand_atoms (
            pdb_file TEXT,
            ligand TEXT,
            chain TEXT,
            sequence_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL,
            virus_name TEXT
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ligand_water_distances (
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
            distance REAL,
            virus_name TEXT
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS solvent_exposed_atoms (
            pdb_file TEXT,
            ligand TEXT,
            chain TEXT,
            sequence_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL,
            virus_name TEXT
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS distal_atoms (
            pdb_file TEXT,
            ligand TEXT,
            chain TEXT,
            sequence_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL,
            virus_name TEXT
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS receptor_binding_pocket (
            pdb_file TEXT,
            residue TEXT,
            chain TEXT,
            atom_id INTEGER,
            exact_atom TEXT,
            atom_type TEXT,
            x REAL,
            y REAL,
            z REAL,
            virus_name TEXT
        );
    ''')

    conn.commit()
    return conn

def insert_csv_data(conn, csv_file, table_name, virus_name):
    cursor = conn.cursor()
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        if reader.fieldnames:  # Ensure file is not empty
            # Rename "PDB File" to "pdb_file" if necessary
            renamed_fieldnames = [col.replace("PDB File", "pdb_file") for col in reader.fieldnames]
            columns = ', '.join([f'[{col}]' for col in renamed_fieldnames] + ["virus_name"])
            placeholders = ', '.join('?' * (len(renamed_fieldnames) + 1))
            query = f"INSERT INTO {table_name} ({columns}) VALUES ({placeholders})"
            
            # Adjust the rows to have "pdb_file" instead of "PDB File"
            rows = [tuple([row.get("PDB File") if "PDB File" in row else row.get("pdb_file")] + list(row.values())[1:] + [virus_name]) for row in reader]
            
            cursor.executemany(query, rows)
    conn.commit()
# Step 3: Load CSV data into the new database format
def load_data_into_database_with_virus_name(base_directory, db_name="viral_data.db"):
    conn = initialize_database_with_virus_name(db_name)
    base_path = Path(base_directory)

    for virus_dir in base_path.iterdir():
        if virus_dir.is_dir():
            virus_name = virus_dir.name  # Get the virus name from the directory
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

    print(f"All data loaded into {db_name} with virus names.")
    conn.close()

# Example usage
base_directory = r"C:\Users\joeys\OneDrive\Documents\Projects\VIRAL_DATABASE\Processed_PDBs"

load_data_into_database_with_virus_name(base_directory)
