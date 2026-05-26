import sqlite3

def extract_pdb_code(conn, table_name):
    cursor = conn.cursor()
    try:
        # Select all rows in the table to process the pdb_file column
        cursor.execute(f"SELECT rowid, pdb_file FROM {table_name}")
        rows = cursor.fetchall()
        
        for row in rows:
            rowid, pdb_file = row
            # Extract the PDB code by splitting the path and removing the extension
            pdb_code = pdb_file.split("\\")[-1].replace(".pdb", "")
            
            # Update the row with the new pdb_code
            cursor.execute(f"UPDATE {table_name} SET pdb_file = ? WHERE rowid = ?", (pdb_code, rowid))

        # Commit changes after processing all rows
        conn.commit()
        print(f"PDB codes updated successfully in {table_name}")
    except Exception as e:
        print(f"Error processing table {table_name}: {e}")

def update_all_tables(db_name="viral_data.db"):
    conn = sqlite3.connect(db_name)
    
    # List of tables to update
    tables = ["distal_atoms", "ligand_atoms", "ligand_water_distances", "solvent_exposed_atoms", "receptor_binding_pocket"]
    
    for table in tables:
        extract_pdb_code(conn, table)
    
    conn.close()
    print("All tables have been updated.")

# Example usage
update_all_tables()
