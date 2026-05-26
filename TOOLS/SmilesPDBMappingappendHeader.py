import sqlite3

# Connect to your SQLite database
conn = sqlite3.connect('viral_data.db')
cursor = conn.cursor()

# Create a new table with the updated headers
cursor.execute('''
    CREATE TABLE IF NOT EXISTS SMILES_MAP_PDB_NEW (
        virus_name TEXT,
        pdb_id TEXT,
        ligand_id TEXT,
        chain TEXT,
        exact_atom TEXT,
        atom_id INTEGER,
        atom_index INTEGER,
        smiles_atom_index INTEGER
    )
''')

# Copy the existing data from the old table to the new table with new headers
cursor.execute('''
    INSERT INTO SMILES_MAP_PDB_NEW (virus_name, pdb_id, ligand_id, chain, exact_atom, atom_id, atom_index, smiles_atom_index)
    SELECT "Virus Name", "PDB ID", "Ligand", "Chain", "PDB Atom Name", "PDB Atom ID", "PDB Atom Index", "SMILES Atom Index"
    FROM SMILES_MAP_PDB
''')

# Drop the old table
cursor.execute('DROP TABLE IF EXISTS SMILES_MAP_PDB')

# Rename the new table to the original name
cursor.execute('ALTER TABLE SMILES_MAP_PDB_NEW RENAME TO SMILES_MAP_PDB')

# Commit changes and close the connection
conn.commit()
conn.close()

print("Headers renamed and data copied successfully.")
