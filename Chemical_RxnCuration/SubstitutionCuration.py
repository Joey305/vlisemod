# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'substitution_reactions_functional_groups.txt'

# Define keywords for Substitution Reactions (SN1 and SN2)
substitution_keywords = [
    'alkyl halide', 'alcohol', 'alkane', 'bromo', 'chloro', 'iodo', 'leaving group'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on substitution keywords
substitution_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in substitution_keywords)
]

# Write the filtered Substitution Reaction SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(substitution_patterns)

print(f"Filtered substitution reaction functional groups saved to {output_file}")
