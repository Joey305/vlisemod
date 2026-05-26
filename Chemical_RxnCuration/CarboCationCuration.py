# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'carbocation_rearrangement_functional_groups.txt'

# Define keywords for Carbocation Rearrangements
carbocation_keywords = [
    'alkyl halide', 'alcohol', 'alkene', 'carbocation', 'haloalkane', 'isopropyl', 'tert-butyl'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on carbocation rearrangement keywords
carbocation_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in carbocation_keywords)
]

# Write the filtered Carbocation Rearrangement SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(carbocation_patterns)

print(f"Filtered carbocation rearrangement functional groups saved to {output_file}")
