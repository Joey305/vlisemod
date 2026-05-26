# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'grignard_functional_groups.txt'

# Define keywords for Grignard Reactions
grignard_keywords = [
    'organohalide', 'carbonyl', 'aldehyde', 'ketone', 'ester', 'bromo', 'chloro', 'iodo'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on Grignard keywords
grignard_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in grignard_keywords)
]

# Write the filtered Grignard SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(grignard_patterns)

print(f"Filtered Grignard functional groups saved to {output_file}")
