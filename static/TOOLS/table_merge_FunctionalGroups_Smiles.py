import sqlite3
import pandas as pd

# Connect to the viral_data.db database
conn = sqlite3.connect('viral_data.db')
cursor = conn.cursor()

# Query Ligand_Arp_Diagram and Functional_GROUPED
query_arp = '''SELECT virus_name, pdb_id, ligand, chain, ligand_id FROM Ligand_Arp_Diagram'''
query_grouped = '''SELECT virus_name, pdb_id, ligand, smiles, functional_groups FROM Functional_GROUPED'''

df_arp = pd.read_sql_query(query_arp, conn)
df_grouped = pd.read_sql_query(query_grouped, conn)

# Merge the two dataframes on virus_name, pdb_id, and ligand
merged_df = pd.merge(df_arp, df_grouped, on=['virus_name', 'pdb_id', 'ligand'], how='inner')

# Create the new table
merged_df.to_sql('Ligand_Atoms_Smiles', conn, if_exists='replace', index=False)

# Commit changes and close the connection
conn.commit()
conn.close()

print("New table 'Ligand_Atoms_Smiles' has been created and saved to viral_data.db")
