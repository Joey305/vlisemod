import sqlite3
import csv
from pathlib import Path
import os

def add_ligand_arp_to_db(db_path, txt_file_path, table_name="Ligand_Arp_Diagram"):
    # Connect to the SQLite database (or create it if it doesn't exist)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create a new table for Ligand ARP data
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS {table_name} (
            virus_name TEXT,
            pdb_id TEXT,
            ligand TEXT,
            chain TEXT,
            residue_id TEXT
        )
    ''')

    # Open the Ligand ARP text file and read the data
    with open(txt_file_path, 'r', encoding='utf-8') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        next(csv_reader)  # Skip the header

        # Insert each row into the database table
        for row in csv_reader:
            cursor.execute(f'''
                INSERT INTO {table_name} (virus_name, pdb_id, ligand, chain, residue_id)
                VALUES (?, ?, ?, ?, ?)
            ''', row)

    # Commit the changes and close the connection
    conn.commit()
    conn.close()

    print(f"Data from {txt_file_path} has been successfully added to {db_path} under the table '{table_name}'.")

# Define the path to the Viral_data.db
db_path = './Viral_data.db'

# Define the path to the combined Ligand ARP Diagram file
# Since the file is in the Sorted_PDB_Files directory, we fetch the most recent one
sorted_pdb_dir = Path('./Database_DATA/Sorted_PDB_Files')
arp_files = sorted_pdb_dir.glob('Ligand_Arp_Diagram*.txt')
latest_arp_file = max(arp_files, key=os.path.getctime)  # Find the latest ARP file by creation time

# Call the function to add the ARP data to the database
add_ligand_arp_to_db(db_path, latest_arp_file)
