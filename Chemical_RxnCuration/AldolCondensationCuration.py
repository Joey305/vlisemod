# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'aldol_condensation_functional_groups.txt'

# Define keywords that are likely to indicate groups involved in aldol condensation
aldol_condensation_keywords = [
    'aldehyde', 'ketone', 'carbonyl', 'enol', 'enolate', 'alpha-hydrogen',
    'alpha carbon', 'conjugated ketone', 'beta-hydroxy'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on aldol condensation keywords
aldol_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in aldol_condensation_keywords)
]

# Write the filtered aldol condensation SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(aldol_patterns)

print(f"Filtered aldol condensation functional groups saved to {output_file}")
