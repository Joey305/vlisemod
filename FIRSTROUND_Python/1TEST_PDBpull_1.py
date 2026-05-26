import requests
from pathlib import Path
import os

def download_pdb_files(virus_name, output_dir):
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
            break

        result_json = response.json()
        pdb_ids = [entry['identifier'] for entry in result_json.get('result_set', [])]
        all_pdb_ids.extend(pdb_ids)
        
        if len(pdb_ids) < rows:
            break  # No more results to retrieve
        
        start += rows

    output_dir_path = Path(output_dir) / virus_name
    output_dir_path.mkdir(parents=True, exist_ok=True)

    for pdb_id in all_pdb_ids:
        pdb_url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        pdb_response = requests.get(pdb_url)
        if pdb_response.status_code == 200:
            with open(output_dir_path / f'{pdb_id}.pdb', 'w') as file:
                file.write(pdb_response.text)
        else:
            print(f"Failed to download {pdb_id}.pdb, status code: {pdb_response.status_code}")

# Example virus list and directory
current_dir = '/Users/jerrismacbook/Downloads/VIRAL_DATABASE/Database_DATA'
virus_list = ['Human immunodeficiency virus 1']

for virus in virus_list:
    download_pdb_files(virus, current_dir)
