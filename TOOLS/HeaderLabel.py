import sqlite3

def update_column_name(db_path, table_name, old_column, new_column):
    # Connect to the database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Check if the old column exists in the table
    cursor.execute(f"PRAGMA table_info({table_name});")
    columns = cursor.fetchall()

    column_names = [col[1] for col in columns]
    if old_column in column_names:
        # Rename the column
        cursor.execute(f"ALTER TABLE {table_name} RENAME COLUMN {old_column} TO {new_column};")
        print(f"Column '{old_column}' has been renamed to '{new_column}' in the '{table_name}' table.")
    else:
        print(f"Column '{old_column}' not found in the '{table_name}' table.")

    # Commit changes and close the connection
    conn.commit()
    conn.close()

# Path to your database
db_path = 'viral_data.db'

# Update the 'sequence_id' column to 'atom_id' in the 'ligand_atoms' table
update_column_name(db_path, 'Functional_GROUPED', 'pdb_file', 'pdb_id')
