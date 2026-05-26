# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'michael_addition_functional_groups.txt'

# Define keywords for Michael Addition
michael_keywords = [
    'α,β-unsaturated carbonyl', 'enolate', 'carbonyl', 'conjugated', 'double bond'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on Michael Addition keywords
michael_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in michael_keywords)
]

# Write the filtered Michael Addition SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(michael_patterns)

print(f"Filtered Michael addition functional groups saved to {output_file}")
