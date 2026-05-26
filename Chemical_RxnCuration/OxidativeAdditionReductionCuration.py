# Define the path to your SMARTS patterns file
input_file = 'functional_groups_smarts.txt'
output_file = 'oxidative_addition_reductive_elimination_functional_groups.txt'

# Define keywords for Oxidative Addition and Reductive Elimination
oxidative_reductive_keywords = [
    'transition metal', 'organometallic', 'metal', 'palladium', 'platinum', 'rhodium', 'nickel'
]

# Read the SMARTS patterns from the input file
with open(input_file, 'r') as file:
    smarts_patterns = file.readlines()

# Filter the SMARTS patterns based on Oxidative Addition and Reductive Elimination keywords
oxidative_reductive_patterns = [
    pattern for pattern in smarts_patterns
    if any(keyword in pattern.lower() for keyword in oxidative_reductive_keywords)
]

# Write the filtered Oxidative Addition and Reductive Elimination SMARTS patterns to the output file
with open(output_file, 'w') as file:
    file.writelines(oxidative_reductive_patterns)

print(f"Filtered oxidative addition/reductive elimination functional groups saved to {output_file}")
