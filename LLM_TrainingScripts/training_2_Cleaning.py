import json
import pandas as pd
import re
from collections import defaultdict

# Input and output file names
INPUT_FILENAME = "training_data2.json"
OUTPUT_FILENAME = "training_data3.json"

# Load your data
with open(INPUT_FILENAME, "r", encoding="utf-8") as f:
    training_data = json.load(f)

# Function to replace unicode escape Angstrom with real symbol
def replace_unicode_angstrom(data):
    if isinstance(data, list):
        return [replace_unicode_angstrom(item) for item in data]
    elif isinstance(data, dict):
        return {k: replace_unicode_angstrom(v) for k, v in data.items()}
    elif isinstance(data, str):
        return data.replace("\u00c5", "Å")
    return data

# Regex to extract interaction components
ligand_pdb_pattern = re.compile(r"What interactions does (.+?) have in PDB (\w{4})\?")
output_pattern = re.compile(r"residue (\w+ \d+) \(chain (\w+)\) via a ([\w_]+) at ([\d.]+) Å")

# Prepare mapping
interaction_map = defaultdict(list)
non_interactions = []

for entry in training_data:
    input_text = entry.get("input", "")
    output_text = entry.get("output", "")

    match = ligand_pdb_pattern.match(input_text)
    if match:
        ligand_raw, pdb_id = match.groups()
        for res, chain, interaction_type, distance in output_pattern.findall(output_text):
            interaction_map[(ligand_raw, pdb_id)].append((res, chain, interaction_type, distance))
    else:
        non_interactions.append(entry)

# Generate merged interaction entries
aggregated_data = []
for (ligand_raw, pdb_id), interactions in interaction_map.items():
    merged = "\n".join([
        f" - Residue {res} (chain {chain}) via {itype} at {dist} Å."
        for res, chain, itype, dist in interactions
    ])
    combined_output = f"Ligand {ligand_raw} in PDB {pdb_id} interacts with:\n{merged}"
    aggregated_data.append({
        "input": f"What interactions does {ligand_raw} have in PDB {pdb_id}?",
        "output": combined_output
    })

# Combine and replace unicode
final_data = replace_unicode_angstrom(non_interactions + aggregated_data)

# Save to new file
with open(OUTPUT_FILENAME, "w", encoding="utf-8") as f:
    json.dump(final_data, f, indent=4, ensure_ascii=False)

print(f"✅ Aggregation complete: saved to {OUTPUT_FILENAME}")
