"""
README: PDB Downloader Script for Single-Strain Virus-Specific PDB Files

DESCRIPTION:
This script is designed to download Protein Data Bank (PDB) files related to viruses with a single strain, such as SARS-CoV-2 (Severe Acute Respiratory Syndrome Coronavirus 2). 
The script queries the RCSB PDB database using the virus name and retrieves all available PDB IDs associated with the virus, downloading the corresponding PDB files. 
The script organizes the PDB files into directories named after each virus.

USAGE SCENARIO:
This script is ideal for downloading PDB files for viruses with a single strain that is heavily studied, such as SARS-CoV-2. 
By focusing on one strain of the virus, researchers can collect all associated PDB files for detailed structural analysis, computational studies, or drug discovery efforts.

INPUT:
- virus_list: A list containing the virus names (strings) for which you want to retrieve PDB files. In this script, we focus on viruses with a single strain.
  Example:
  virus_list = [
      "Severe acute respiratory syndrome coronavirus 2"
  ]

- output_dir: The base directory where all PDB files will be stored. For each virus in the `virus_list`, a subdirectory will be created within this base directory, provided that PDB files are found.
  Example:
  output_dir = 'Database_DATA'

OUTPUT:
- The script will generate a subdirectory for each virus that has associated PDB files. These subdirectories will be created within the specified `output_dir`.
  For example, if the virus "Severe acute respiratory syndrome coronavirus 2" has associated PDB files, a directory named:
  `Database_DATA/Severe acute respiratory syndrome coronavirus 2/`
  will be created.

- Inside the virus-specific subdirectory, the corresponding PDB files will be saved. For instance, the subdirectory might contain files such as:
  - '6M17.pdb'
  - '7K3F.pdb'

  The files are named using their respective PDB ID, which is standard practice for naming PDB files.

HOW TO USE:
1. **Modify the `virus_list` variable**: 
   - Add the names of the viruses you are studying, particularly if they have a single dominant strain. 
   - In the current version, the script is set up for "Severe acute respiratory syndrome coronavirus 2" (SARS-CoV-2), but you can replace this with any other virus name that has a unique strain.

2. **Set the `output_dir` variable**: 
   - Specify the base directory where you want to save the downloaded PDB files.
   - Example: `output_dir = 'Database_DATA'`

3. **Run the script**: 
   - The script will search for all PDB files associated with the virus names provided in the `virus_list` and download them to the specified output directory. Subdirectories will be created dynamically as needed.

OUTPUT USAGE:
- **Research and Analysis**: The downloaded PDB files can be used for various types of structural analysis, such as:
  - Molecular docking studies
  - Molecular dynamics simulations
  - Protein-ligand interactions
  - Drug discovery and design
  - Protein structure visualization and annotation

- **Focus on Single-Strain Viruses**: This script is particularly useful for viruses that have a single dominant strain under study, like SARS-CoV-2. It allows researchers to gather all relevant structural data for a specific strain in one organized directory.

KEY VARIABLES:
- `virus_list`: A list of virus names (strings) for which PDB files should be downloaded. Customize this list to include the virus names you are researching.
  Example: `["Severe acute respiratory syndrome coronavirus 2"]`

- `output_dir`: The directory where the downloaded PDB files will be saved. The script will create subdirectories for each virus within this directory.
  Example: `output_dir = 'Database_DATA'`

DEPENDENCIES:
- Python 3.x
- `requests` library (can be installed with `pip install requests`)
- `pathlib` and `os` are part of the standard Python library.

NOTES:
- **Single-Strain Focus**: The script is tailored to viruses with a single strain that are studied intensively. However, it can also be modified to handle multi-strain viruses if needed by adjusting the virus naming in `virus_list`.
- **Progress Reporting**: The script provides progress updates every 60 seconds while downloading PDB files. This helps users monitor the download progress, especially for large datasets.
- **Automatic Folder Creation**: A folder will only be created if PDB files are found for the virus. Empty folders are not generated, ensuring that only relevant data is stored.
"""

import requests
from pathlib import Path
import os
import time
import logging


current_dir = 'Database_DATA'
virus_list = [
    "Severe acute respiratory syndrome coronavirus 2"
    # "Human immunodeficiency virus 1"

]









def setup_logging(virus_name, output_dir):
    log_dir = Path(output_dir) / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f'{virus_name}_download.log'
    
    logging.basicConfig(filename=log_file, 
                        filemode='w', 
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    
    logging.info(f"Starting download process for {virus_name}")
    
def download_pdb_files(virus_name, output_dir):
    setup_logging(virus_name, output_dir)
    
    all_pdb_ids = []
    start = 0
    rows = 100
    while True:
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                    "operator": "exact_match",
                    "value": virus_name
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": start,
                    "rows": rows
                },
                "scoring_strategy": "combined"
            }
        }

        url = 'https://search.rcsb.org/rcsbsearch/v2/query'
        response = requests.post(url, json=query)
        
        # Debugging output to understand the response
        print(f"Response status code: {response.status_code}")
        try:
            result_json = response.json()
            print(f"Response JSON: {result_json}")
        except ValueError:
            print("Failed to decode JSON response")
            logging.error("Failed to decode JSON response")
            break
        
        if response.status_code != 200:
            print(f"Failed to retrieve PDB IDs for {virus_name}, status code: {response.status_code}")
            logging.error(f"Failed to retrieve PDB IDs for {virus_name}, status code: {response.status_code}")
            break

        pdb_ids = [entry['identifier'] for entry in result_json.get('result_set', [])]
        all_pdb_ids.extend(pdb_ids)
        
        if len(pdb_ids) < rows:
            break  # No more results to retrieve
        
        start += rows

    total_pdb_files = len(all_pdb_ids)
    print(f"Total number of PDB files found: {total_pdb_files}")
    logging.info(f"Total number of PDB files found: {total_pdb_files}")

    output_dir_path = Path(output_dir) / virus_name
    output_dir_path.mkdir(parents=True, exist_ok=True)

    downloaded_count = 0
    skipped_count = 0
    start_time = time.time()

    for pdb_id in all_pdb_ids:
        pdb_file_path = output_dir_path / f'{pdb_id}.pdb'
        
        # Check if the file already exists
        if pdb_file_path.exists():
            print(f"Skipping {pdb_id}.pdb, already exists.")
            logging.info(f"Skipping {pdb_id}.pdb, already exists.")
            skipped_count += 1
            continue

        pdb_url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        pdb_response = requests.get(pdb_url)
        if pdb_response.status_code == 200:
            with open(pdb_file_path, 'w') as file:
                file.write(pdb_response.text)
            downloaded_count += 1
            logging.info(f"Successfully downloaded {pdb_id}.pdb")
        else:
            print(f"Failed to download {pdb_id}.pdb, status code: {pdb_response.status_code}")
            logging.error(f"Failed to download {pdb_id}.pdb, status code: {pdb_response.status_code}")

        # Print update every 60 seconds
        if time.time() - start_time > 60:
            print(f"Downloaded {downloaded_count} of {total_pdb_files} PDB files (Skipped {skipped_count})")
            logging.info(f"Downloaded {downloaded_count} of {total_pdb_files} PDB files (Skipped {skipped_count})")
            start_time = time.time()

    # Final update after completion
    print(f"Finished downloading {downloaded_count} of {total_pdb_files} PDB files (Skipped {skipped_count})")
    logging.info(f"Finished downloading {downloaded_count} of {total_pdb_files} PDB files (Skipped {skipped_count})")


for virus in virus_list:
    print(f"Trying search term: {virus}")
    download_pdb_files(virus, current_dir)
