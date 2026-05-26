import sqlite3

# Connect to your SQLite database
conn = sqlite3.connect('viral_data.db')
cursor = conn.cursor()

# Rename the RUPLEY_SASA_DATA table with the updated headers
cursor.execute('''
    CREATE TABLE IF NOT EXISTS RUPLEY_SASA_DATA_NEW (
        virus_name TEXT,
        pdb_id TEXT,
        ligand_id TEXT,
        chain TEXT,
        exact_atom TEXT,
        atom_id INTEGER,
        SASA_Area REAL
    )
''')

# Copy the existing data from the old table to the new table with new headers
cursor.execute('''
    INSERT INTO RUPLEY_SASA_DATA_NEW (virus_name, pdb_id, ligand_id, chain, exact_atom, atom_id, SASA_Area)
    SELECT "Virus Name", "PDB ID", "Ligand", "Chain", "Atom Name", "Atom ID", "SASA Value"
    FROM RUPLEY_SASA_DATA
''')

# Drop the old table
cursor.execute('DROP TABLE IF EXISTS RUPLEY_SASA_DATA')

# Rename the new table to the original name
cursor.execute('ALTER TABLE RUPLEY_SASA_DATA_NEW RENAME TO RUPLEY_SASA_DATA')

# Commit changes and close the connection
conn.commit()
conn.close()

print("Headers renamed and data copied successfully.")
