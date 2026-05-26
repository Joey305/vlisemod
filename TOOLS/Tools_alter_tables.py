import sqlite3

def add_virus_name_column(db_name="viral_data.db"):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # List of tables to alter
    tables = ["ligand_atoms", "ligand_water_distances", "solvent_exposed_atoms", "distal_atoms", "receptor_binding_pocket"]

    for table in tables:
        try:
            cursor.execute(f"ALTER TABLE {table} ADD COLUMN virus_name TEXT;")
            print(f"Added 'virus_name' column to {table}")
        except sqlite3.OperationalError as e:
            # This might happen if the column already exists
            print(f"Error altering table {table}: {e}")

    conn.commit()
    conn.close()

# Run the script
add_virus_name_column()
