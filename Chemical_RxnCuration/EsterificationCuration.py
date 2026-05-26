# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'esterification_functional_groups.txt'

# Define keywords for Esterification
esterification_keywords = [
    'carboxylic acid', 'alcohol', 'acid chloride', 'anhydride', 'hydroxyl', 'ester linkage'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on esterification keywords
esterification_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in esterification_keywords)
]

# Write the filtered esterification SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(esterification_patterns)

print(f"Filtered esterification functional groups saved to {output_file}")
