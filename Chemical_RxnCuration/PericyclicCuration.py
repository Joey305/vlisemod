# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'pericyclic_functional_groups.txt'

# Define keywords for Pericyclic Reactions
pericyclic_keywords = [
    'conjugated diene', 'alkyne', 'alkene', 'cyclic', 'double bond', 'diels-alder'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on pericyclic keywords
pericyclic_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in pericyclic_keywords)
]

# Write the filtered Pericyclic SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(pericyclic_patterns)

print(f"Filtered pericyclic functional groups saved to {output_file}")
