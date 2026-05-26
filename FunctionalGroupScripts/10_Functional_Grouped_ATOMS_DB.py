import os
import csv
from pathlib import Path

# Function to merge all files into one with the additional virus name column
def merge_functional_group_atoms_files(base_input_dir, virus_names, output_file):
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        # Write the header to the output file
        header = ['virus_name', 'PDB ID', 'Ligand ID', 'Chain', 'Functional Group', 'Atom ID', 'Exact Atom', 'Atom Type']
        writer.writerow(header)

        for virus_name in virus_names:
            input_file = Path(base_input_dir) / 'Sorted_PDB_Files' / virus_name / 'with_ligands_smiles_functional_groups-reorganized.txt'

            if not input_file.exists():
                print(f"File not found for virus '{virus_name}': {input_file}")
                continue

            with open(input_file, 'r') as infile:
                reader = csv.reader(infile, delimiter='\t')
                next(reader)  # Skip the header row

                for row in reader:
                    # Prepend the virus name to the row
                    new_row = [virus_name] + row
                    writer.writerow(new_row)

            print(f"Merged data from virus '{virus_name}' into the output file.")

    print(f"All files have been merged into {output_file}")

# Main function
def main():
    base_input_dir = './Database_DATA'
    output_file = 'merged_functional_groups_atoms.txt'

    virus_names = [
        'Human immunodeficiency virus 1',
        "Severe acute respiratory syndrome coronavirus 2",
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
        "Human papillomavirus type 6a"
    ]

    # Merge files
    merge_functional_group_atoms_files(base_input_dir, virus_names, output_file)

if __name__ == "__main__":
    main()
