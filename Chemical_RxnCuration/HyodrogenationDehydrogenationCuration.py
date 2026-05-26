# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'hydrogenation_dehydrogenation_functional_groups.txt'

# Define keywords for Hydrogenation and Dehydrogenation
hydrogenation_keywords = [
    'alkene', 'alkyne', 'double bond', 'carbonyl', 'reduction', 'hydride', 'catalyst', 'hydrogen'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on Hydrogenation and Dehydrogenation keywords
hydrogenation_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in hydrogenation_keywords)
]

# Write the filtered Hydrogenation and Dehydrogenation SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(hydrogenation_patterns)

print(f"Filtered hydrogenation/dehydrogenation functional groups saved to {output_file}")
