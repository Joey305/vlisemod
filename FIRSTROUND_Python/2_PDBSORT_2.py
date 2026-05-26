"""
# PDB Parsing Script for Viral Data with Ligand Detection

## Overview
This script is designed to parse PDB files for a specified virus and classify them based on the presence of ligand codes. It identifies whether a PDB file contains ligands of interest and separates them into two categories: `with_ligands.csv` for files containing ligands of interest and `without_ligands.csv` for those containing only non-interest ligands (defined by a preset list).

The script dynamically names the output CSV files using the current date and time, ensuring no overwrites occur when the script is executed multiple times.

## Purpose
- **Ligand Detection**: Identify PDB files with ligands of interest (i.e., ligands not included in a preset exclusion list).
- **Classification**: Separate PDB files into two CSV files: `with_ligands.csv` for files containing ligands of interest and `without_ligands.csv` for those containing only non-interest ligands.
- **Automated Output Naming**: The output CSV files are named with the current timestamp to prevent overwriting older outputs.

## Prerequisites
To run the script, you must have:
1. **Python 3.x**: Ensure Python is installed on your system.
2. **CSV and PDB Files**: The script expects PDB files to be organized in folders specific to a virus type, and it processes only the files within that virus's directory.
3. **Non-ligand Codes List**: The script includes a preset list of ligand codes that are treated as non-interest ligands.

## Directory Structure
The script expects the following directory structure:

./Database_DATA └── Virus_Name/ └── PDB_Files (e.g., *.pdb)


Where `Virus_Name` is the name of the virus, and `PDB_Files` are the PDB files associated with that virus. For example:
./Database_DATA/Human immunodeficiency virus 1/1A8G.pdb



## Output
The output CSV files will be saved in the following directory structure:
./Database_DATA/Sorted_PDB_Files/Virus_Name/ └── with_ligands_YYYYMMDD_HHMMSS.csv └── without_ligands_YYYYMMDD_HHMMSS.csv


Where `YYYYMMDD_HHMMSS` is the current date and time at which the script was executed.

### Output File Contents
1. **with_ligands.csv**
   - This file contains the PDB IDs and ligand codes for files that include ligands of interest.
   - Format: `PDB ID, Ligand Codes`
   
2. **without_ligands.csv**
   - This file contains the PDB IDs and any non-interest ligands (e.g., water, ions).
   - Format: `PDB ID, Ligand Codes (Non-interest)`

## Script Parameters
### Non-ligand Codes
The script includes a list of non-interest ligands (e.g., water, ions, common solvents), which are excluded from the list of ligands of interest. If a PDB file contains only these non-interest ligands, it will be categorized into the `without_ligands.csv` file.

### Virus Name
The variable `virus_name` should be set to the name of the virus you want to process. This is used to locate the appropriate input folder and store the output files in the correct directory.

### Base Input Directory
`base_input_dir` is the base directory where all virus-specific folders are stored. This defaults to:
./Database_DATA/



### Base Output Directory
`base_output_dir` is the base directory where the script stores the resulting CSV files. This defaults to:
./Database_DATA/Sorted_PDB_Files/



## How to Use the Script

### 1. Prepare the Input Data
Ensure that your input PDB files are organized in a directory structure like:
./Database_DATA/<Virus_Name>/*.pdb


For example:
./Database_DATA/Human immunodeficiency virus 1/1A8G.pdb



### 2. Set the Virus Name
Update the `virus_name` variable in the script to match the name of the virus you are processing. For example:
```python
virus_name = 'Human immunodeficiency virus 1'
3. Run the Script
From the directory containing the script, execute the following command:

python <script_name>.py
This will parse the PDB files in the specified virus folder, generate the necessary CSV files, and store them in the appropriate output directory.

4. Check the Output
After the script runs, you will find two CSV files in the output directory:

with_ligands_<timestamp>.csv: Contains the PDB files with ligands of interest.
without_ligands_<timestamp>.csv: Contains the PDB files that only have non-interest ligands.
Notes
If no ligands of interest are found, the script will still generate the without_ligands.csv file with the relevant non-interest ligands.
Ensure that the PDB files are correctly formatted for proper ligand detection. """



import os
import csv
from pathlib import Path
from datetime import datetime

def parse_pdb_files(virus_name, base_input_dir, base_output_dir):
    non_ligand_codes = set([
        "03N", "1GP", "1PE", "12P", "2CO", "2DA", "2FX", "2GL", "2DX", "3DR", "3PA", "41H", "4BA", "43X", "4DX", "73O", "90A", "ABA", "ACA", "ACE", "ACY", "ACT", "ADN", 
        "ADP", "ADE", "AIB", "AMP", "ANP", "ARF", "ARS", "ASN", "ATM", "AZT", "ATP", "B3E", "B3P", "BAM", "BDG", "BDF", "BEN", "BGC", "BIF", "BMA", "BME", 
        "BOG", "BR", "BTE", "BU3", "BYZ", "CAC", "CAF", "CAS", "CA", "CD", "CIT", "CL", "CME", "CMO", "CMM", "CS", "CSD", "CSO", "CO", "CO3", "CU", 
        "DAL", "DAS", "DBU", "DDG", "DHI", "DIV", "DIQ", "DLE", "DLY", "DMS", "DCY", "DGL", "DGN", "DPN", "DOC", "DOD", "DIL", "DPV", "DTD", "DTR", 
        "DTP", "DPP", "DTT", "DTY", "DVA", "DPR", "DTV", "EDO", "EPE", "ESD", "FAD", "FLC", "FE", "FG7", "FMT", "FRU", "FUC", "FUM", "G46", "G47", "GGL", 
        "GAL", "G3p", "GLO", "GLC", "GLY", "GOA", "GVE", "GOL", "GLU", "HAI", "HEM", "HEP", "HEZ", "HG", "HOH", "HPH", "IDG", "IIL", "IMD", "IOD", "IPA", 
        "IVA", "K", "KCX", "KDO", "KF2", "MAN", "LAC", "LEU", "LYS", "MES", "MEA", "MG", "MLA", "MN", "M7G", "MK8", "MNK", "MPD", "MPT", "MRG", "MSE", "MYR", "NAD", "NA", 
        "NAG", "NAO", "NEN", "NH2", "NH4", "NI", "NLE", "NHE", "NO3", "NTB", "OAS", "OIL", "OMC", "OXY", "PB", "PC", "PCA", "PEG", "PGE", "PG2", "PG3", "PG4", 
        "P6G", "PH2", "PUT", "PO4", "P03", "P6S", "PPI", "PTR", "PYZ", "QNC", "RIB", "SIA", "SLZ", "SNU", "SMC", "SO3", "SO2",  "SO4", "STA", "SRT", "STA", "TAR", "TAM", 
        "TEO", "TPO", "TRS", "TYS", "U2X", "UB4", "U3X",  "UZ1", "UZ4", "UZ7", "UMP", "URE", "UNX", "VLM", "VME", "XCP", "XPC", "YCM", "ZN", "Z9N"
    ])

    # Create output directory path for each virus
    output_dir_path = Path(base_output_dir) / virus_name
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Get the current date to dynamically name the output files
    current_date = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Define paths for output CSV and TXT files
    with_ligands_path = output_dir_path / f'with_ligands_{current_date}.csv'
    without_ligands_path = output_dir_path / f'without_ligands_{current_date}.csv'
    with_ligands_txt_path = output_dir_path / f'with_ligands_{current_date}.txt'
    without_ligands_txt_path = output_dir_path / f'without_ligands_{current_date}.txt'
    with_ligands_arp_txt_path = output_dir_path / f'with_ligands_ARP_{current_date}.txt'

    # Initialize CSV writers and TXT writers
    with open(with_ligands_path, 'w', newline='', encoding='utf-8') as with_ligands_file, \
         open(without_ligands_path, 'w', newline='', encoding='utf-8') as without_ligands_file, \
         open(with_ligands_txt_path, 'w', encoding='utf-8') as with_ligands_txt_file, \
         open(without_ligands_txt_path, 'w', encoding='utf-8') as without_ligands_txt_file, \
         open(with_ligands_arp_txt_path, 'w', encoding='utf-8') as with_ligands_arp_txt_file:
        
        # Write headers
        with_ligands_writer = csv.writer(with_ligands_file)
        without_ligands_writer = csv.writer(without_ligands_file)
        with_ligands_writer.writerow(['PDB ID', 'Ligand Codes'])
        without_ligands_writer.writerow(['PDB ID', 'Ligand Codes (Non-interest)'])
        
        with_ligands_txt_file.write("PDB ID\tLigand Codes\n")
        without_ligands_txt_file.write("PDB ID\tLigand Codes (Non-interest)\n")
        with_ligands_arp_txt_file.write("PDB ID\tLigand\tChain\tResidue_ID\n")
        
        # Iterate over PDB files in the input directory
        input_dir_path = Path(base_input_dir) / virus_name
        for pdb_file in input_dir_path.glob('*.pdb'):
            pdb_id = pdb_file.stem
            has_ligand_of_interest = False
            ligand_codes = set()
            non_ligand_found = set()
            unique_ligand_info = set()  # Track unique ligand-chain-residue combinations
            
            with pdb_file.open('r') as file:
                for line in file:
                    if line.startswith('HETATM'):
                        hetatm_code = line[17:20].strip()
                        chain_id = line[21].strip()  # Chain ID
                        residue_id = line[22:26].strip()  # Residue ID
                        
                        if hetatm_code not in non_ligand_codes:
                            has_ligand_of_interest = True
                            ligand_codes.add(hetatm_code)
                            unique_ligand_info.add((pdb_id, hetatm_code, chain_id, residue_id))
                        else:
                            non_ligand_found.add(hetatm_code)
            
            # Convert ligand codes and non-ligand codes to strings
            ligand_codes_str = ', '.join(ligand_codes)
            non_ligand_found_str = ', '.join(non_ligand_found)

            # Write the data to the CSV and TXT files
            if has_ligand_of_interest:
                with_ligands_writer.writerow([pdb_id, ligand_codes_str])
                with_ligands_txt_file.write(f"{pdb_id}\t{ligand_codes_str}\n")
                
                # Write the unique ligand-chain-residue information once for each unique combination
                for info in unique_ligand_info:
                    with_ligands_arp_txt_file.write(f"{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}\n")
            else:
                without_ligands_writer.writerow([pdb_id, non_ligand_found_str if non_ligand_found else 'None'])
                without_ligands_txt_file.write(f"{pdb_id}\t{non_ligand_found_str if non_ligand_found else 'None'}\n")

# Define base directories relative to the current script location
base_input_dir = './Database_DATA'
base_output_dir = './Database_DATA/Sorted_PDB_Files'

# List of viruses to process
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

# Run the parsing function for each virus in the list
for virus in virus_list:
    print(f"Processing virus: {virus}")
    parse_pdb_files(virus, base_input_dir, base_output_dir)




# import os
# import csv
# from pathlib import Path
# from datetime import datetime

# def parse_pdb_files(virus_name, base_input_dir, base_output_dir):
#     non_ligand_codes = set([
#         "03N", "1GP", "1PE", "12P", "2CO", "2DA", "2FX", "2GL", "2DX", "3DR", "3PA", "41H", "4BA", "43X", "4DX", "73O", "90A", "ABA", "ACA", "ACE", "ACY", "ACT", "ADN", 
#         "ADP", "ADE", "AIB", "AMP", "ANP", "ARF", "ARS", "ASN", "ATM", "AZT", "ATP", "B3E", "B3P", "BAM", "BDG", "BDF", "BEN", "BGC", "BIF", "BMA", "BME", 
#         "BOG", "BR", "BTE", "BU3", "BYZ", "CAC", "CAF", "CAS", "CA", "CD", "CIT", "CL", "CME", "CMO", "CMM", "CS", "CSD", "CSO", "CO", "CO3", "CU", 
#         "DAL", "DAS", "DBU", "DDG", "DHI", "DIV", "DIQ", "DLE", "DLY", "DMS", "DCY", "DGL", "DGN", "DPN", "DOC", "DOD", "DIL", "DPV", "DTD", "DTR", 
#         "DTP", "DPP", "DTT", "DTY", "DVA", "DPR", "DTV", "EDO", "EPE", "ESD", "FAD", "FLC", "FE", "FG7", "FMT", "FRU", "FUC", "FUM", "G46", "G47", "GGL", 
#         "GAL", "G3p", "GLO", "GLC", "GLY", "GOA", "GVE", "GOL", "GLU", "HAI", "HEM", "HEP", "HEZ", "HG", "HOH", "HPH", "IDG", "IIL", "IMD", "IOD", "IPA", 
#         "IVA", "K", "KCX", "KDO", "KF2", "MAN", "LAC", "LEU", "LYS", "MES", "MEA", "MG", "MLA", "MN", "M7G", "MK8", "MNK", "MPD", "MPT", "MRG", "MSE", "MYR", "NAD", "NA", 
#         "NAG", "NAO", "NEN", "NH2", "NH4", "NI", "NLE", "NHE", "NO3", "NTB", "OAS", "OIL", "OMC", "OXY", "PB", "PC", "PCA", "PEG", "PGE", "PG2", "PG3", "PG4", 
#         "P6G", "PH2", "PUT", "PO4", "P03", "P6S", "PPI", "PTR", "PYZ", "QNC", "RIB", "SIA", "SLZ", "SNU", "SMC", "SO3", "SO2",  "SO4", "STA", "SRT", "STA", "TAR", "TAM", 
#         "TEO", "TPO", "TRS", "TYS", "U2X", "UB4", "U3X",  "UZ1", "UZ4", "UZ7", "UMP", "URE", "UNX", "VLM", "VME", "XCP", "XPC", "YCM", "ZN", "Z9N"
#     ])
#     # Create output directory path for each virus
#     output_dir_path = Path(base_output_dir) / virus_name
#     output_dir_path.mkdir(parents=True, exist_ok=True)
    
#     # Get the current date to dynamically name the CSV files
#     current_date = datetime.now().strftime("%Y%m%d_%H%M%S")
    
#     # Define paths for output CSV files with dynamic names
#     with_ligands_path = output_dir_path / f'with_ligands_{current_date}.csv'
#     without_ligands_path = output_dir_path / f'without_ligands_{current_date}.csv'
    
#     # Initialize CSV writers
#     with open(with_ligands_path, 'w', newline='', encoding='utf-8') as with_ligands_file, \
#          open(without_ligands_path, 'w', newline='', encoding='utf-8') as without_ligands_file:
        
#         with_ligands_writer = csv.writer(with_ligands_file)
#         without_ligands_writer = csv.writer(without_ligands_file)
        
#         # Write headers
#         with_ligands_writer.writerow(['PDB ID', 'Ligand Codes'])
#         without_ligands_writer.writerow(['PDB ID', 'Ligand Codes (Non-interest)'])
        
#         # Iterate over PDB files in the input directory
#         input_dir_path = Path(base_input_dir) / virus_name
#         for pdb_file in input_dir_path.glob('*.pdb'):
#             pdb_id = pdb_file.stem
#             has_ligand_of_interest = False
#             ligand_codes = set()
#             non_ligand_found = set()
            
#             with pdb_file.open('r') as file:
#                 for line in file:
#                     if line.startswith('HETATM'):
#                         # Ensure that the ligand code is explicitly treated as a string
#                         hetatm_code = line[17:20].strip()
                        
#                         # Check if the ligand code is a non-interest ligand
#                         if hetatm_code not in non_ligand_codes:
#                             has_ligand_of_interest = True
#                             ligand_codes.add(hetatm_code)
#                         else:
#                             non_ligand_found.add(hetatm_code)
            
#             # Convert ligand codes and non-ligand codes to strings
#             ligand_codes_str = ', '.join(ligand_codes)
#             non_ligand_found_str = ', '.join(non_ligand_found)

#             # Write the data to the CSV
#             if has_ligand_of_interest:
#                 with_ligands_writer.writerow([pdb_id, ligand_codes_str])
#             else:
#                 without_ligands_writer.writerow([pdb_id, non_ligand_found_str if non_ligand_found else 'None'])

#     # Define paths for output text files with dynamic names
#     with_ligands_path_txt = output_dir_path / f'with_ligands_{current_date}.txt'
#     without_ligands_path_txt = output_dir_path / f'without_ligands_{current_date}.txt'
    
#     # Open text files for writing
#     with open(with_ligands_path_txt, 'w', encoding='utf-8') as with_ligands_file, \
#          open(without_ligands_path_txt, 'w', encoding='utf-8') as without_ligands_file:
        
#         # Write headers
#         with_ligands_file.write("PDB ID\tLigand Codes\n")
#         without_ligands_file.write("PDB ID\tLigand Codes (Non-interest)\n")
        
#         # Iterate over PDB files in the input directory
#         for pdb_file in input_dir_path.glob('*.pdb'):
#             pdb_id = pdb_file.stem
#             has_ligand_of_interest = False
#             ligand_codes = set()
#             non_ligand_found = set()
            
#             with pdb_file.open('r') as file:
#                 for line in file:
#                     if line.startswith('HETATM'):
#                         hetatm_code = line[17:20].strip()
                        
#                         # Check if the ligand code is a non-interest ligand
#                         if hetatm_code not in non_ligand_codes:
#                             has_ligand_of_interest = True
#                             ligand_codes.add(hetatm_code)
#                         else:
#                             non_ligand_found.add(hetatm_code)
            
#             # Convert ligand codes and non-ligand codes to strings
#             ligand_codes_str = ', '.join(ligand_codes)
#             non_ligand_found_str = ', '.join(non_ligand_found)

#             # Write the data to the text files
#             if has_ligand_of_interest:
#                 with_ligands_file.write(f"{pdb_id}\t{ligand_codes_str}\n")
#             else:
#                 without_ligands_file.write(f"{pdb_id}\t{non_ligand_found_str if non_ligand_found else 'None'}\n")

# # Define base directories relative to the current script location
# base_input_dir = './Database_DATA' # + "/HPV"  #########COMMENTED OUT HPV FOR WHEN WORKING WITH OTHER VIRUSES
# base_output_dir = './Database_DATA/Sorted_PDB_Files'

# # List of viruses to process
# virus_list = [
    

# ######  iF WORKING WITH HPV FILES, MAKE SURE TO 
# ####### COMMENT OUT THE BELOW FILES THAT ARE  
# ####### NOT IN THE DIRECTORY HPV

#     # 'Severe acute respiratory syndrome coronavirus 2',
#     "Human immunodeficiency virus 1"


    
#     #########################################
#     ######IF WORKING WITH HPV, MAKE SURE TO CHANGE 
#     # THE BASE INPUT TO INCLUDE /HPV AT THE END SO 
#     # THAT IT CAN FIND THE BELOW DIRECTORIES THAT 
#     # WERE SAVED UNDER THE BLANKET NAME OF HPV TO 
#     # HELP WITH ORGANIZATION
#      #############################################
#     # "Human papillomavirus",
#     # "Human papillomavirus 1",
#     # "Human papillomavirus 11",
#     # "Human papillomavirus 16",
#     # "Human papillomavirus 18",
#     # "Human papillomavirus 26",
#     # "Human papillomavirus 31",
#     # "Human papillomavirus 33",
#     # "Human papillomavirus 35",
#     # "Human papillomavirus 4",
#     # "Human papillomavirus 45",
#     # "Human papillomavirus 49",
#     # "Human papillomavirus 51",
#     # "Human papillomavirus 52",
#     # "Human papillomavirus 53",
#     # "Human papillomavirus 58",
#     # "Human papillomavirus 59",
#     # "Human papillomavirus 6",
#     # "Human papillomavirus 66",
#     # "Human papillomavirus type 16",
#     # "Human papillomavirus type 6",
#     # "Human papillomavirus type 6a"



# ]


# # Run the parsing function for each virus in the list
# for virus in virus_list:
#     print(f"Processing virus: {virus}")
#     parse_pdb_files(virus, base_input_dir, base_output_dir)

