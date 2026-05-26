import pandas as pd
import sqlite3

def load_data_to_database(txt_path, db_path, table_name):
    print("Loading data from text file...")
    # Read the data from the text file
    data = []
    with open(txt_path, 'r') as file:
        for line in file:
            parts = line.strip().split(':')
            ligand = parts[0].strip()
            synonym = parts[1].strip() if len(parts) > 1 else ""  # Handle blanks in synonyms
            data.append({'ligand': ligand, 'synonym': synonym})

    # Create a DataFrame from the data
    df = pd.DataFrame(data)
    print("Data loaded into DataFrame with columns:", df.columns)

    # Connect to the SQLite database
    print("Connecting to the database...")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Ensure the table exists, create it if not (assuming basic structure: ligand and synonym)
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS {table_name} (
            ligand TEXT,
            synonym TEXT,
            UNIQUE(ligand, synonym)  -- Ensure unique combinations of ligand and synonym
        )
    ''')

    # Loop through the DataFrame and insert data while checking for duplicates
    for index, row in df.iterrows():
        ligand = row['ligand']
        synonym = row['synonym']

        # Check if the ligand and synonym already exist in the database
        cursor.execute(f'''
            SELECT COUNT(*) FROM {table_name}
            WHERE ligand = ? AND synonym = ?
        ''', (ligand, synonym))
        result = cursor.fetchone()

        # If the combination does not exist, insert the new data
        if result[0] == 0:
            cursor.execute(f'''
                INSERT INTO {table_name} (ligand, synonym)
                VALUES (?, ?)
            ''', (ligand, synonym))
            print(f"Inserted new record: {ligand} - {synonym if synonym else 'No synonym'}")
        else:
            print(f"Record already exists: {ligand} - {synonym if synonym else 'No synonym'}")

    # Commit the transaction and close the database connection
    conn.commit()
    conn.close()
    print("Database connection closed. Data insertion complete.")

# Path to your text file and SQLite database
text_file_path = 'Ligand_Synonyms_Corrected.txt'
database_path = 'viral_data.db'
table_name = 'ligand_synonyms'

# Run the function
load_data_to_database(text_file_path, database_path, table_name)
