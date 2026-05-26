# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'friedel_crafts_functional_groups.txt'

# Define keywords for Friedel-Crafts Alkylation/Acylation
friedel_crafts_keywords = [
    'aromatic', 'benzene', 'toluene', 'phenyl', 'alkyl halide', 'acyl halide', 'chloro',
    'bromo', 'iodo', 'aryl', 'alkyl chloride', 'alkyl bromide', 'alkyl iodide'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on Friedel-Crafts keywords
friedel_crafts_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in friedel_crafts_keywords)
]

# Write the filtered Friedel-Crafts SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(friedel_crafts_patterns)

print(f"Filtered Friedel-Crafts functional groups saved to {output_file}")
