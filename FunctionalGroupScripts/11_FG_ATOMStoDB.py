import sqlite3
import csv
from pathlib import Path

# Function to create a table if it doesn't exist
def create_table_if_not_exists(conn):
    create_table_query = '''
    CREATE TABLE IF NOT EXISTS Functional_Group_Atoms (
        virus_name TEXT,
        pdb_id TEXT,
        ligand_id TEXT,
        chain TEXT,
        functional_group TEXT,
        atom_id TEXT,
        exact_atom TEXT,
        atom_type TEXT
    )
    '''
    conn.execute(create_table_query)
    conn.commit()

# Function to insert data into the database
def insert_data_into_db(conn, csv_file):
    with open(csv_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        next(reader)  # Skip the header row
        
        insert_query = '''
        INSERT INTO Functional_Group_Atoms (virus_name, pdb_id, ligand_id, chain, functional_group, atom_id, exact_atom, atom_type)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        '''
        
        for row in reader:
            conn.execute(insert_query, row)
        
        conn.commit()

# Main function to handle the database connection and data insertion
def main():
    db_path = 'viral_data.db'
    csv_file = 'merged_functional_groups_atoms.txt'

    # Connect to the database
    conn = sqlite3.connect(db_path)

    # Create the table if it doesn't exist
    create_table_if_not_exists(conn)

    # Insert data from the CSV into the database
    insert_data_into_db(conn, csv_file)

    # Close the database connection
    conn.close()
    print(f"Data from {csv_file} has been successfully added to the Functional_Group_Atoms table in the database {db_path}.")

if __name__ == "__main__":
    main()
