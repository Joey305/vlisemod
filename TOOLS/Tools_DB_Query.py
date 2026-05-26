import sqlite3

def query_database(db_name="viral_data.db"):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # Query each table
    tables = ["ligand_atoms", "ligand_water_distances", "solvent_exposed_atoms", "distal_atoms", "receptor_binding_pocket", "Functional_GROUPED", 'Functional_Group_Atoms']

    for table in tables:
        print(f"Contents of {table}:")
        try:
            cursor.execute(f"SELECT * FROM {table} LIMIT 10")  # Fetch first 10 rows to check the data
            rows = cursor.fetchall()
            if rows:
                for row in rows:
                    print(row)
            else:
                print(f"No data found in {table}")
        except sqlite3.OperationalError as e:
            print(f"Error querying table {table}: {e}")
        print("\n")

    conn.close()

# Example usage
query_database("viral_data.db")
