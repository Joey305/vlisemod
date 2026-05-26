"""
README: CIF Downloader Script for Virus-Specific PDB Files
DESCRIPTION:
This script allows you to download Protein Data Bank (PDB) files in CIF format for a given virus or set of viruses by querying the RCSB PDB database. It retrieves all available PDB IDs associated with the specified viruses and downloads the corresponding CIF files.

Originally used for handling large numbers of Human Papillomavirus (HPV) strains, the script has since been modified to accommodate a wider range of viruses. The downloaded CIF files are stored in a single directory for simplicity, regardless of the number of viruses or strains.

INPUT:
virus_list: A list of virus names (strings). The script queries the PDB database for each virus in the list and attempts to download all associated CIF files. Example:

python
Copy code
virus_list = [
    "Human papillomavirus",
    "Human immunodeficiency virus 1"
]
output_dir: The base directory where the downloaded CIF files will be stored. All CIF files for the viruses in the list will be saved in this directory. Example: output_dir = 'Database_DATA/All_Viruses'

OUTPUT:
All CIF files for the viruses in the virus_list are stored together in the specified output_dir. The script no longer creates subdirectories for each virus but instead stores all files in the same directory.

Example structure:

Copy code
Database_DATA/All_Viruses/
  ├── 1ABC.cif
  ├── 2DEF.cif
  └── 3XYZ.cif
HOW TO USE:
Customize the virus_list variable to include the names of the viruses you want to search for in the PDB database. You can add as many virus names as necessary.
Set the output_dir variable to the desired directory where the CIF files will be saved.
Run the script. It will query the RCSB PDB database, download all CIF files related to the viruses in the virus_list, and store them in the specified output_dir.
KEY VARIABLES:
virus_list: List of virus names to query. Modify this list with the specific viruses you are interested in. The script will search the RCSB PDB database for exact matches of the names in this list. Example:

python
Copy code
virus_list = ["Human papillomavirus", "Human immunodeficiency virus 1"]
output_dir: The base directory where downloaded CIF files will be stored. All virus-related CIF files will be stored together in this directory. Example: "Database_DATA/All_Viruses"

DEPENDENCIES:
Python 3.x
requests library (can be installed with pip install requests)
pathlib and os are part of the standard Python library.
NOTES:
The script only downloads CIF files if they exist for a given virus. If no CIF files are found, no files will be saved.
The script prints status updates every 60 seconds during the download process to keep the user informed of the progress.
"""


import requests
from pathlib import Path
import os
import time
import logging

# Step 1: Setup logging for the download process
def setup_logging(virus_name, output_dir):
    log_dir = Path(output_dir) / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f'{virus_name}_download.log'
    
    logging.basicConfig(filename=log_file, 
                        filemode='w', 
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    
    logging.info(f"Starting download process for {virus_name}")

# Step 2: Download CIF files from RCSB PDB
def download_cif_files(virus_name, output_dir):
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

        if response.status_code != 200:
            print(f"Failed to retrieve PDB IDs for {virus_name}, status code: {response.status_code}")
            logging.error(f"Failed to retrieve PDB IDs for {virus_name}, status code: {response.status_code}")
            break

        try:
            result_json = response.json()
        except ValueError:
            print("Failed to decode JSON response")
            logging.error("Failed to decode JSON response")
            break

        pdb_ids = [entry['identifier'] for entry in result_json.get('result_set', [])]
        all_pdb_ids.extend(pdb_ids)
        
        if len(pdb_ids) < rows:
            break  # No more results to retrieve
        
        start += rows

    total_pdb_files = len(all_pdb_ids)
    print(f"Total number of CIF files found: {total_pdb_files}")
    logging.info(f"Total number of CIF files found: {total_pdb_files}")

    if total_pdb_files == 0:
        print(f"No CIF files found for {virus_name}. Skipping...")
        logging.info(f"No CIF files found for {virus_name}. Skipping...")
        return

    # Set up the output directory for Arpeggio_Data
    output_dir_path = Path(output_dir)  / virus_name
    output_dir_path.mkdir(parents=True, exist_ok=True)

    downloaded_count = 0
    start_time = time.time()

    for pdb_id in all_pdb_ids:
        cif_url = f'https://files.rcsb.org/download/{pdb_id}.cif'
        cif_file_path = output_dir_path / f'{pdb_id}.cif'

        # Download the CIF file
        cif_response = requests.get(cif_url)
        if cif_response.status_code == 200:
            with open(cif_file_path, 'w') as file:
                file.write(cif_response.text)
            logging.info(f"Successfully downloaded {pdb_id}.cif")
            downloaded_count += 1
        else:
            print(f"Failed to download {pdb_id}.cif, status code: {cif_response.status_code}")
            logging.error(f"Failed to download {pdb_id}.cif, status code: {cif_response.status_code}")

        # Print update every 60 seconds
        if time.time() - start_time > 60:
            print(f"Downloaded {downloaded_count} of {total_pdb_files} CIF files")
            logging.info(f"Downloaded {downloaded_count} of {total_pdb_files} CIF files")
            start_time = time.time()

    # Final update after completion
    print(f"Finished downloading {downloaded_count} of {total_pdb_files} CIF files")
    logging.info(f"Finished downloading {downloaded_count} of {total_pdb_files} CIF files")

# Step 3: Specify the directory and virus list
current_dir = 'Arpeggio_DATA/'
virus_list = virus_list = [
    "Severe acute respiratory syndrome coronavirus 2",
    "Human immunodeficiency virus 1",
    "HPV",
    "Human papillomavirus",
    "Human papillomavirus 1",
    "Human papillomavirus 2",
    "Human papillomavirus 3",
    "Human papillomavirus 4",
    "Human papillomavirus 5",
    "Human papillomavirus 6",
    "Human papillomavirus 7",
    "Human papillomavirus 8",
    "Human papillomavirus 9",
    "Human papillomavirus 10",
    "Human papillomavirus 11",
    "Human papillomavirus 12",
    "Human papillomavirus 13",
    "Human papillomavirus 14",
    "Human papillomavirus 15",
    "Human papillomavirus 16",
    "Human papillomavirus 17",
    "Human papillomavirus 18",
    "Human papillomavirus 19",
    "Human papillomavirus 20",
    "Human papillomavirus 21",
    "Human papillomavirus 22",
    "Human papillomavirus 23",
    "Human papillomavirus 24",
    "Human papillomavirus 25",
    "Human papillomavirus 26",
    "Human papillomavirus 27",
    "Human papillomavirus 28",
    "Human papillomavirus 29",
    "Human papillomavirus 30",
    "Human papillomavirus 31",
    "Human papillomavirus 32",
    "Human papillomavirus 33",
    "Human papillomavirus 34",
    "Human papillomavirus 35",
    "Human papillomavirus 36",
    "Human papillomavirus 37",
    "Human papillomavirus 38",
    "Human papillomavirus 39",
    "Human papillomavirus 40",
    "Human papillomavirus 41",
    "Human papillomavirus 42",
    "Human papillomavirus 43",
    "Human papillomavirus 44",
    "Human papillomavirus 45",
    "Human papillomavirus 46",
    "Human papillomavirus 47",
    "Human papillomavirus 48",
    "Human papillomavirus 49",
    "Human papillomavirus 50",
    "Human papillomavirus 51",
    "Human papillomavirus 52",
    "Human papillomavirus 53",
    "Human papillomavirus 54",
    "Human papillomavirus 55",
    "Human papillomavirus 56",
    "Human papillomavirus 57",
    "Human papillomavirus 58",
    "Human papillomavirus 59",
    "Human papillomavirus 60",
    "Human papillomavirus 61",
    "Human papillomavirus 62",
    "Human papillomavirus 63",
    "Human papillomavirus 64",
    "Human papillomavirus 65",
    "Human papillomavirus 66",
    "Human papillomavirus 67",
    "Human papillomavirus 68",
    "Human papillomavirus 69",
    "Human papillomavirus 70",
    "Human papillomavirus 71",
    "Human papillomavirus 72",
    "Human papillomavirus 73",
    "Human papillomavirus 74",
    "Human papillomavirus 75",
    "Human papillomavirus 76",
    "Human papillomavirus 77",
    "Human papillomavirus 78",
    "Human papillomavirus 79",
    "Human papillomavirus 80",
    "Human papillomavirus 81",
    "Human papillomavirus 82",
    "Human papillomavirus 83",
    "Human papillomavirus 84",
    "Human papillomavirus 85",
    "Human papillomavirus 86",
    "Human papillomavirus 87",
    "Human papillomavirus 88",
    "Human papillomavirus 89",
    "Human papillomavirus 90",
    "Human papillomavirus 91",
    "Human papillomavirus 92",
    "Human papillomavirus 93",
    "Human papillomavirus 94",
    "Human papillomavirus 95",
    "Human papillomavirus 96",
    "Human papillomavirus 97",
    "Human papillomavirus 98",
    "Human papillomavirus 99",
    "Human papillomavirus 100",
    "Human papillomavirus 101",
    "Human papillomavirus 102",
    "Human papillomavirus 103",
    "Human papillomavirus 104",
    "Human papillomavirus 105",
    "Human papillomavirus 106",
    "Human papillomavirus 107",
    "Human papillomavirus 108",
    "Human papillomavirus 109",
    "Human papillomavirus 110",
    "Human papillomavirus 111",
    "Human papillomavirus 112",
    "Human papillomavirus 113",
    "Human papillomavirus 114",
    "Human papillomavirus 115",
    "Human papillomavirus 116",
    "Human papillomavirus 117",
    "Human papillomavirus 118",
    "Human papillomavirus 119",
    "Human papillomavirus 120",
    "Human papillomavirus 121",
    "Human papillomavirus 122",
    "Human papillomavirus 123",
    "Human papillomavirus 124",
    "Human papillomavirus 125",
    "Human papillomavirus 126",
    "Human papillomavirus 127",
    "Human papillomavirus 128",
    "Human papillomavirus 129",
    "Human papillomavirus 130",
    "Human papillomavirus 131",
    "Human papillomavirus 132",
    "Human papillomavirus 133",
    "Human papillomavirus 134",
    "Human papillomavirus 135",
    "Human papillomavirus 136",
    "Human papillomavirus 137",
    "Human papillomavirus 138",
    "Human papillomavirus 139",
    "Human papillomavirus 140",
    "Human papillomavirus 141",
    "Human papillomavirus 142",
    "Human papillomavirus 143",
    "Human papillomavirus 144",
    "Human papillomavirus 145",
    "Human papillomavirus 146",
    "Human papillomavirus 147",
    "Human papillomavirus 148",
    "Human papillomavirus 149",
    "Human papillomavirus 150",
    "Human papillomavirus 151",
    "Human papillomavirus 152",
    "Human papillomavirus 153",
    "Human papillomavirus 154",
    "Human papillomavirus 155",
    "Human papillomavirus 156",
    "Human papillomavirus 157",
    "Human papillomavirus 158",
    "Human papillomavirus 159",
    "Human papillomavirus 160",
    "Human papillomavirus type 6",
    "Human papillomavirus type 6a",
    "Human papillomavirus type 16",
]


# Step 4: Download CIF files for each virus in the list
for virus in virus_list:
    print(f"Trying search term: {virus}")
    download_cif_files(virus, current_dir)
