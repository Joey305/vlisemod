import sqlite3
import pandas as pd

# Define the database name
db_name = "viral_data.db"

def csv_to_sqlite(csv_file, table_name, db_name):
    """
    Loads a CSV file into an SQLite table.
    
    Parameters:
    csv_file (str): Path to the CSV file.
    table_name (str): Name of the SQLite table to store the data.
    db_name (str): Name of the SQLite database.
    """
    # Read the CSV into a pandas DataFrame
    df = pd.read_csv(csv_file)
    
    # Connect to the SQLite database (or create it if it doesn't exist)
    conn = sqlite3.connect(db_name)
    
    # Write the data to the specified SQLite table
    df.to_sql(table_name, conn, if_exists='replace', index=False)
    
    # Commit the changes and close the connection
    conn.commit()
    conn.close()
    
    print(f"Data from {csv_file} successfully written to {table_name} table in {db_name}.")

# File paths for the two CSV files
rupley_sasa_csv = "Rupley_SASA_data_batch.csv"
pdb_smile_code_csv = "PDB_SMILE_Code_Batch.csv"

# Define the table names for the database
rupley_sasa_table = "RUPLEY_SASA_DATA"
smiles_map_pdb_table = "SMILES_MAP_PDB"

# Write the CSV files to the SQLite database
csv_to_sqlite(rupley_sasa_csv, rupley_sasa_table, db_name)
csv_to_sqlite(pdb_smile_code_csv, smiles_map_pdb_table, db_name)
