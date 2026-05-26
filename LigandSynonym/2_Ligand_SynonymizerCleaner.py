import pandas as pd

def clean_ligand_synonyms(input_file, output_file):
    # Read the original Ligand_Synonyms.txt file
    df = pd.read_csv(input_file, sep='\t')

    # Replace "No synonyms data available" with an empty string, keeping the ligand intact
    df['Synonyms'] = df['Synonyms'].replace("No synonyms data available", "")

    # Save the corrected data to a new file
    df.to_csv(output_file, index=False, sep='\t')
    print(f"Corrected synonyms text file has been created successfully as {output_file}.")

# Paths to your input and output files
input_file_path = 'Ligand_Synonyms.txt'
output_file_path = 'Ligand_Synonyms_Corrected.txt'

# Call the function to clean the file
clean_ligand_synonyms(input_file_path, output_file_path)
