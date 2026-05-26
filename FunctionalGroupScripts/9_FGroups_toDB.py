import sqlite3
import csv
from pathlib import Path

# Function to create a table if it doesn't exist
def create_table_if_not_exists(conn):
    create_table_query = '''
    CREATE TABLE IF NOT EXISTS Functional_GROUPED (
        virus_name TEXT,
        pdb_file TEXT,
        ligand TEXT,
        smiles TEXT,
        functional_groups TEXT
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
        INSERT INTO Functional_GROUPED (virus_name, pdb_file, ligand, smiles, functional_groups)
        VALUES (?, ?, ?, ?, ?)
        '''
        
        for row in reader:
            conn.execute(insert_query, row)
        
        conn.commit()

# Main function to handle the database connection and data insertion
def main():
    db_path = 'viral_data.db'
    csv_file = 'merged_grouped_functional_groups.txt'

    # Connect to the database
    conn = sqlite3.connect(db_path)

    # Create the table if it doesn't exist
    create_table_if_not_exists(conn)

    # Insert data from the CSV into the database
    insert_data_into_db(conn, csv_file)

    # Close the database connection
    conn.close()
    print(f"Data from {csv_file} has been successfully added to the Functional_GROUPED table in the database {db_path}.")

if __name__ == "__main__":
    main()
