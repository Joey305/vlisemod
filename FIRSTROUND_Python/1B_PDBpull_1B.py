"""
README: PDB Downloader Script for Virus-Specific PDB Files

DESCRIPTION:
This script allows you to download Protein Data Bank (PDB) files for a given virus or set of viruses by querying the RCSB PDB database. 
It retrieves all available PDB IDs associated with the virus and downloads the corresponding PDB files. 
The script organizes the PDB files into directories named after each virus for which PDB files were found.

********************************************************************************************************************************
********************************************************************************************************************************
*******THIS SCRIPT WAS ALTERED TO SPECIFICALLY ACCOMODATE THE LARGE NUMBER OF VIRUS STRAINS FOR HPV*****************************
*******PLEASE ONLY USE THIS VERSION OF THE SCRIPT WHEN WORKING WITH A SINGLE VIRUS THAT CONTAINS MULTIPLE***********************
*******STRAINS. OTHER WISE YOU CAN PULL AND ACCOMODATE MULTIPLE DIFFERENT TYPES OF VIRUSES WITH THE 1A_PDBpull_.py**************
********************************************************************************************************************************
********************************************************************************************************************************
INPUT:
- virus_list: A list of virus names (strings). The script will query the PDB database for each virus in the list and attempt to download all associated PDB files.
  Example:
  virus_list = [
      "Human papillomavirus",
      "Human immunodeficiency virus 1"
  ]

- output_dir: The base directory where the downloaded PDB files will be stored. The script will create a subdirectory for each virus within this base directory, but only if there are PDB files to download.
  Example:
  output_dir = 'Database_DATA/HPV'

OUTPUT:
- For each virus that has associated PDB files, the script will create a folder named after the virus inside the specified output directory. 
  Example folder structure (if PDB files exist for these viruses):
  Database_DATA/HPV/Human papillomavirus/
  Database_DATA/HPV/Human papillomavirus 1/

  Inside each virus folder, the corresponding PDB files will be stored. 
  Example file inside a folder: '1ABC.pdb', '2DEF.pdb'

HOW TO USE:
1. Customize the `virus_list` variable to include the names of the viruses you want to search for in the PDB database. You can add as many virus names as necessary.
2. Set the `output_dir` variable to the desired base directory where the PDB files will be saved.
3. Run the script. It will query the RCSB PDB database, download all PDB files related to the viruses in the `virus_list`, and organize them by virus name.

KEY VARIABLES:
- `virus_list`: List of virus names to query. Modify this list with the specific viruses you are interested in. The script will search the RCSB PDB database for exact matches of the names in this list.
  Example: ["Human papillomavirus", "Human immunodeficiency virus 1"]
  
- `output_dir`: The base directory where downloaded PDB files will be stored. Subdirectories for each virus with PDB files will be created within this directory.
  Example: "Database_DATA/HPV"

DEPENDENCIES:
- Python 3.x
- `requests` library (can be installed with `pip install requests`)
- `pathlib` and `os` are part of the standard Python library.

NOTES:
- The script only creates subdirectories for viruses that actually have PDB files available for download. If no PDB files are found for a virus, no directory will be created.
- The script prints status updates every 60 seconds during downloading to keep the user informed of the progress.

"""


import requests
from pathlib import Path
import os
import time
import logging

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

    if total_pdb_files == 0:
        print(f"No PDB files found for {virus_name}. Skipping...")
        logging.info(f"No PDB files found for {virus_name}. Skipping...")
        return

    output_dir_path = Path(output_dir) / virus_name
    output_dir_path.mkdir(parents=True, exist_ok=True)

    downloaded_count = 0
    start_time = time.time()

    for pdb_id in all_pdb_ids:
        pdb_url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        pdb_response = requests.get(pdb_url)
        if pdb_response.status_code == 200:
            with open(output_dir_path / f'{pdb_id}.pdb', 'w') as file:
                file.write(pdb_response.text)
            downloaded_count += 1
            logging.info(f"Successfully downloaded {pdb_id}.pdb")
        else:
            print(f"Failed to download {pdb_id}.pdb, status code: {pdb_response.status_code}")
            logging.error(f"Failed to download {pdb_id}.pdb, status code: {pdb_response.status_code}")

        # Print update every 60 seconds
        if time.time() - start_time > 60:
            print(f"Downloaded {downloaded_count} of {total_pdb_files} PDB files")
            logging.info(f"Downloaded {downloaded_count} of {total_pdb_files} PDB files")
            start_time = time.time()

    # Final update after completion
    print(f"Finished downloading {downloaded_count} of {total_pdb_files} PDB files")
    logging.info(f"Finished downloading {downloaded_count} of {total_pdb_files} PDB files")


current_dir = 'Database_Data/'
virus_list = [
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

for virus in virus_list:
    print(f"Trying search term: {virus}")
    download_pdb_files(virus, current_dir)
   
