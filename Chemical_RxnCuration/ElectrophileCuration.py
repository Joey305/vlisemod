# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'electrophilic_addition_functional_groups.txt'

# Define keywords that are likely to indicate electrophilic addition groups
electrophilic_addition_keywords = [
    'alkene', 'alkyne', 'arene', 'aromatic', 'double bond', 'triple bond', 
    'olefin', 'vinyl', 'benzene', 'phenyl', 'conjugated', 'diene', 'aromatic ring'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on electrophilic addition keywords
electrophilic_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in electrophilic_addition_keywords)
]

# Write the filtered electrophilic addition SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(electrophilic_patterns)

print(f"Filtered electrophilic addition functional groups saved to {output_file}")
