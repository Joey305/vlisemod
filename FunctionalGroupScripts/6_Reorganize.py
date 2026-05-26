import csv
from collections import defaultdict

def reorganize_functional_groups(input_file, output_file):
    # Read the input file and categorize data by functional group
    functional_groups = defaultdict(list)
    
    with open(input_file, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        
        for row in reader:
            key = (row['PDB ID'], row['Ligand ID'], row['Chain'], row['Functional Group'])
            functional_groups[key].append(row)

    # Sort the appended rows back into the appropriate functional group
    sorted_functional_groups = {}
    
    for key, rows in functional_groups.items():
        sorted_functional_groups[key] = sorted(rows, key=lambda x: int(x['Atom ID']))
    
    # Write the sorted data to the output file
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        
        for key in sorted_functional_groups:
            for row in sorted_functional_groups[key]:
                writer.writerow(row)

def main():
    base_input_dir = './Database_DATA'
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
        # Add more virus names as needed
    ]

    for virus_name in virus_names:
        input_file = f"{base_input_dir}/Sorted_PDB_Files/{virus_name}/with_ligands_smiles_functional_groups-appended.txt"
        output_file = input_file.replace('-appended.txt', '-reorganized.txt')
        
        reorganize_functional_groups(input_file, output_file)
        print(f"Processed and reorganized functional groups for virus: {virus_name}")

if __name__ == "__main__":
    main()


