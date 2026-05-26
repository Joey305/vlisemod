import os
from pathlib import Path
from datetime import datetime

def append_ligand_arp_files(virus_list, base_output_dir, final_output_file):
    # Get the current date to dynamically name the final output file
    current_date = datetime.now().strftime("%Y%m%d_%H%M%S")
    final_output_path = Path(base_output_dir) / f'Ligand_Arp_Diagram_{current_date}.txt'
    
    # Open the final output file in append mode
    with open(final_output_path, 'w', encoding='utf-8') as final_output:
        final_output.write("Virus Name\tPDB ID\tLigand\tChain\tResidue_ID\n")  # Write the header
        
        # Loop through each virus
        for virus_name in virus_list:
            virus_dir = Path(base_output_dir) / virus_name
            
            # Check if the virus directory exists
            if virus_dir.exists() and virus_dir.is_dir():
                # Look for files that match 'with_ligands_ARP_*.txt'
                for arp_file in virus_dir.glob('with_ligands_ARP_*.txt'):
                    with open(arp_file, 'r', encoding='utf-8') as f:
                        # Skip the header of the ARP file (assumes the first line is the header)
                        next(f)
                        # Read each line, prepend the virus name, and write to the final output file
                        for line in f:
                            final_output.write(f"{virus_name}\t{line}")
            else:
                print(f"Directory for virus '{virus_name}' not found. Skipping...")

    print(f"All ligand ARP files have been appended to {final_output_path}.")

# Define base directory where the sorted PDB files are located
base_output_dir = './Database_DATA/Sorted_PDB_Files'

# List of viruses to process (same as the virus list from the previous script)
virus_list = [
    "Human immunodeficiency virus 1",
    "Human papillomavirus",
    "Human papillomavirus 1",
    "Human papillomavirus 11",
    "Human papillomavirus 16",
    "Human papillomavirus 18",
    "Human papillomavirus 26",
    "Human papillomavirus 31",
    "Human papillomavirus 33",
    "Human papillomavirus 35",
    "Human papillomavirus 4",
    "Human papillomavirus 45",
    "Human papillomavirus 49",
    "Human papillomavirus 51",
    "Human papillomavirus 52",
    "Human papillomavirus 53",
    "Human papillomavirus 58",
    "Human papillomavirus 59",
    "Human papillomavirus 6",
    "Human papillomavirus 66",
    "Human papillomavirus type 16",
    "Human papillomavirus type 6",
    "Human papillomavirus type 6a",
    'Severe acute respiratory syndrome coronavirus 2',
]

# Call the function to append all ARP files
append_ligand_arp_files(virus_list, base_output_dir, 'Ligand_Arp_Diagram.txt')
