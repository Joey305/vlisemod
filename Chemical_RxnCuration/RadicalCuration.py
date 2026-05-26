# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'radical_reaction_functional_groups.txt'

# Define keywords that are likely to indicate radical reaction groups
radical_reaction_keywords = [
    'alkyl halide', 'halogen', 'benzyl', 'allyl', 'alkene', 'aromatic', 'peroxide', 
    'arene', 'double bond', 'alkyne', 'nitro', 'phenyl', 'conjugated', 'diene', 'benzylic'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on radical reaction keywords
radical_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in radical_reaction_keywords)
]

# Write the filtered radical reaction SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(radical_patterns)

print(f"Filtered radical reaction functional groups saved to {output_file}")
