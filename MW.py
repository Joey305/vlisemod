###ADD MOLECULAR WEIGHT OF VIRAL LIGANDS TO DATABASE TO STORE FOR DOWNSTREAM USAGE IN molecular_weight 

import sqlite3
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Use the correct database file
db_path = "viral_data.db"  # Ensure this matches your actual database name
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Check if the table exists
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = [t[0] for t in cursor.fetchall()]
print("Tables in the database:", tables)

if "Ligand_Atoms_Smiles" not in tables:
    print("❌ Error: The table 'Ligand_Atoms_Smiles' does not exist in viral_data.db.")
    conn.close()
    exit()

# Load Ligand_Atoms_Smiles table
query = 'SELECT * FROM "Ligand_Atoms_Smiles"'
ligand_df = pd.read_sql_query(query, conn)

# Check if 'molecular_weight' column exists; if not, add it
cursor.execute("PRAGMA table_info('Ligand_Atoms_Smiles')")
columns = [col[1] for col in cursor.fetchall()]
if "molecular_weight" not in columns:
    cursor.execute("ALTER TABLE 'Ligand_Atoms_Smiles' ADD COLUMN molecular_weight REAL")
    conn.commit()
    print("✅ Added 'molecular_weight' column to Ligand_Atoms_Smiles table.")

# Function to calculate molecular weight from SMILES
def calculate_molecular_weight(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return round(Descriptors.MolWt(mol), 2)  # Rounded for consistency
    except Exception as e:
        print(f"⚠️ Error calculating molecular weight for SMILES {smiles}: {e}")
    return None

# Apply function to compute molecular weight
ligand_df["molecular_weight"] = ligand_df["smiles"].apply(calculate_molecular_weight)

# Update database with the new molecular weight values
for index, row in ligand_df.iterrows():
    cursor.execute(
        "UPDATE 'Ligand_Atoms_Smiles' SET molecular_weight = ? WHERE virus_name = ? AND pdb_id = ? AND ligand = ? AND chain = ? AND ligand_id = ?",
        (row["molecular_weight"], row["virus_name"], row["pdb_id"], row["ligand"], row["chain"], row["ligand_id"]),
    )

# Commit changes and close connection
conn.commit()
conn.close()
print("✅ Molecular weight values successfully added to the database!")
