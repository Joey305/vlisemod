import sqlite3
import requests
from bs4 import BeautifulSoup
import pandas as pd

def fetch_and_clean_ligand_synonyms(db_path):
    print("Connecting to database...")
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Query to select unique ligands
    query = "SELECT DISTINCT ligand FROM Ligand_Arp_Diagram"
    cursor.execute(query)
    ligands = cursor.fetchall()
    print(f"Found {len(ligands)} unique ligands.")
    
    # Close the database connection
    conn.close()
    
    # Prepare to collect data
    ligand_data = []
    clean_rows = []

    # Iterate over each ligand to fetch its synonyms from RCSB PDB
    for (ligand,) in ligands:
        print(f"Fetching data for ligand: {ligand}")
        url = f"https://www.rcsb.org/ligand/{ligand}"
        response = requests.get(url)
        
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            synonym_section = soup.find('tr', id='chemicalSynonyms')
            
            if synonym_section:
                synonyms = synonym_section.find('td').text.strip() if synonym_section.find('td') else "No synonyms found"
                print(f"Synonyms found for {ligand}: {synonyms}")
            else:
                synonyms = "No synonyms data available"
                print(f"No synonyms section found for {ligand}.")
        else:
            synonyms = "Failed to fetch data"
            print(f"Error: Could not retrieve data from RCSB PDB for {ligand}. Status code: {response.status_code}")
        
        ligand_data.append({'Ligand': ligand, 'Synonyms': synonyms})
        
        # Process synonyms for clean data
        if synonyms != "No synonyms data available" and synonyms != "Failed to fetch data":
            synonym_list = [syn.strip() for syn in synonyms.split(';') if syn.strip()]
            for synonym in synonym_list:
                clean_rows.append(f"{ligand}: {synonym}")

    # Save the raw data to a text file
    df = pd.DataFrame(ligand_data)
    df.to_csv('Ligand_Synonyms.txt', index=False, sep='\t')
    print("Raw synonyms text file has been created successfully.")
    
    # Save the cleaned data to a new text file
    with open('Ligand_Synonyms_Clean.txt', 'w') as file:
        for row in clean_rows:
            file.write(row + '\n')
    print("Cleaned synonyms text file has been created successfully.")

# Path to your SQLite database
database_path = 'viral_data.db'
fetch_and_clean_ligand_synonyms(database_path)
