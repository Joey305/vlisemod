# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'hydrolysis_functional_groups.txt'

# Define keywords that indicate groups involved in hydrolysis reactions
hydrolysis_keywords = [
    'ester', 'amide', 'anhydride', 'carbamate', 'lactone', 'lactam', 'phosphate ester',
    'sulfate ester', 'phosphoric acid', 'carboxylic acid', 'ester linkage', 'carbonyl'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on hydrolysis keywords
hydrolysis_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in hydrolysis_keywords)
]

# Write the filtered hydrolysis SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(hydrolysis_patterns)

print(f"Filtered hydrolysis functional groups saved to {output_file}")
