# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'carbanion_reaction_functional_groups.txt'

# Define keywords for Carbanion Reactions (Enolate Chemistry)
carbanion_keywords = [
    'carbonyl', 'aldehyde', 'ketone', 'ester', 'enolate', 'carbonyl group'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on Carbanion keywords
carbanion_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in carbanion_keywords)
]

# Write the filtered Carbanion Reaction SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(carbanion_patterns)

print(f"Filtered carbanion reaction functional groups saved to {output_file}")
