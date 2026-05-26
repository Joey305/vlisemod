# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'diels_alder_functional_groups.txt'

# Define keywords that indicate groups involved in Diels-Alder reactions
diels_alder_keywords = [
    'conjugated diene', 'diene', 'dienophile', 'alkene', 'alkyne', 
    'electron withdrawing', 'carbonyl', 'cyano', 'ester', 'nitro'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on Diels-Alder keywords
diels_alder_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in diels_alder_keywords)
]

# Write the filtered Diels-Alder SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(diels_alder_patterns)

print(f"Filtered Diels-Alder functional groups saved to {output_file}")
