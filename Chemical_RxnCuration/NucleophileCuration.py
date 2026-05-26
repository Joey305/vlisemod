# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'nucleophilic_functional_groups.txt'

# Define keywords that are likely to indicate nucleophilic groups
nucleophilic_keywords = [
    'amine', 'thiol', 'sulfide', 'alcohol', 'alkoxide', 'cyanide', 'imide', 
    'amide', 'sulfhydryl', 'enolate', 'azide', 'nitrogen', 'sulfur', 'sulfate', 
    'hydroxyl', 'methoxy', 'carbamate'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on nucleophilic keywords
nucleophilic_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in nucleophilic_keywords)
]

# Write the filtered nucleophilic SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(nucleophilic_patterns)

print(f"Filtered nucleophilic functional groups saved to {output_file}")
