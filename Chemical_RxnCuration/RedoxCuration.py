# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'redox_functional_groups.txt'

# Define keywords that indicate groups involved in redox reactions
redox_keywords = [
    'alcohol', 'aldehyde', 'ketone', 'carboxylic acid', 'quinone', 
    'thiol', 'amine', 'sulfide', 'disulfide', 'nitro', 'metal',
    'peroxide', 'hydroperoxide', 'epoxide', 'ether', 'phosphine',
    'sulfate', 'phosphate', 'carbonyl', 'cyanide', 'reduction', 'oxidation'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on redox keywords
redox_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in redox_keywords)
]

# Write the filtered redox SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(redox_patterns)

print(f"Filtered redox functional groups saved to {output_file}")
