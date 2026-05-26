# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'decarboxylation_functional_groups.txt'

# Define keywords for Decarboxylation
decarboxylation_keywords = [
    'carboxylic acid', 'carboxylate', 'carbon dioxide', 'beta-keto acid', 'malonic acid'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on decarboxylation keywords
decarboxylation_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in decarboxylation_keywords)
]

# Write the filtered decarboxylation SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(decarboxylation_patterns)

print(f"Filtered decarboxylation functional groups saved to {output_file}")
