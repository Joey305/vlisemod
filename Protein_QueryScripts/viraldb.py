import sqlite3
import csv

# Define the database and CSV file paths
db_file = 'viral_data.db'
csv_file = 'protein_summaryCIF.csv'

# Connect to the database
conn = sqlite3.connect(db_file)
cursor = conn.cursor()

# Create the Virus_Proteins table if it doesn't exist
cursor.execute('''
CREATE TABLE IF NOT EXISTS Virus_Proteins (
    virus_name TEXT,
    pdb_id TEXT,
    protein TEXT
)
''')

# Open the CSV file and insert its content into the Virus_Proteins table
with open(csv_file, 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        cursor.execute('''
        INSERT INTO Virus_Proteins (virus_name, pdb_id, protein)
        VALUES (?, ?, ?)
        ''', (row['virus_name'], row['pdb_id'], row['protein']))

# Commit the transaction and close the connection
conn.commit()
conn.close()

print("Data from protein_summaryCIF.csv has been successfully added to the Virus_Proteins table.")
