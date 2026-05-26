import sqlite3
import csv
import os

# Database connection
db_path = "viral_data.db"
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# CSV file path
csv_file = "Updated_Ligand_Interaction_DiagramData.csv"

# Check if the Arpeggio_Contacts_Data table exists, if not create it with the additional Atom_ID column
cursor.execute(''' 
    CREATE TABLE IF NOT EXISTS Arpeggio_Contacts_Data (
        Virus_Name TEXT,
        PDB_ID TEXT,
        Ligand TEXT,
        Ligand_Number TEXT,
        Chain TEXT,
        Contact TEXT,
        Distance REAL,
        Ligand_Atom TEXT,
        Atom_ID TEXT,  -- New column for Atom_ID
        Residue TEXT,
        Residue_Number INTEGER,
        Residue_Atom TEXT,
        Residue_Chain TEXT
    )
''')
conn.commit()

# Insert data from the CSV file into the Arpeggio_Contacts_Data table
def insert_csv_into_database(csv_file, cursor):
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            cursor.execute('''
                INSERT INTO Arpeggio_Contacts_Data (
                    Virus_Name, PDB_ID, Ligand, Ligand_Number, Chain, Contact, Distance, 
                    Ligand_Atom, Atom_ID, Residue, Residue_Number, Residue_Atom, Residue_Chain
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                row['Virus_Name'], row['PDB_ID'], row['Ligand'], row['Ligand_Number'], row['Chain'],
                row['Contact'], row['Distance'], row['Ligand_Atom'], row['Atom_ID'],  # Including Atom_ID here
                row['Residue'], row['Residue_Number'], row['Residue_Atom'], row['Residue_Chain']
            ))
        conn.commit()

# Insert the data from CSV into the database
insert_csv_into_database(csv_file, cursor)

# Close the database connection
conn.close()

print(f"Data from {csv_file} has been successfully inserted into the Arpeggio_Contacts_Data table in {db_path}.")
