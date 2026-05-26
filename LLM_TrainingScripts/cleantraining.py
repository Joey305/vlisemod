# import json
# import pandas as pd

# def replace_unicode_angstrom(data):
#     """Recursively replaces Unicode escape for Angstrom (\u00c5) with the actual symbol in all string values."""
#     if isinstance(data, list):
#         return [replace_unicode_angstrom(item) for item in data]
#     elif isinstance(data, dict):
#         return {k: replace_unicode_angstrom(v) for k, v in data.items()}
#     elif isinstance(data, str):
#         return data.replace("\u00c5", "Å")
#     return data

# # Load existing training data
# with open("training_data.json", "r", encoding="utf-8") as f:
#     training_data = json.load(f)

# # Convert training data into a DataFrame for easier manipulation
# df = pd.DataFrame(training_data)

# # Separate interaction-related data
# interaction_df = df[df["input"].str.startswith("What interactions does")].copy()

# # Extract ligand and PDB details separately
# interaction_df[["ligand"]] = interaction_df["input"].str.extract(r"What interactions does (\w+) have")
# interaction_df[["pdb_id"]] = interaction_df["input"].str.extract(r"PDB (\w+)\?")

# # Extract residue and chain
# interaction_df[["residue", "chain"]] = interaction_df["output"].str.extract(r"residue (\w+ \d+) \(chain (\w+)\)")

# # Extract interaction type and distance
# interaction_df[["interaction_type", "distance"]] = interaction_df["output"].str.extract(r"via a (\w+) at ([\d.]+) Å\.")

# # Drop rows where necessary data is missing
# interaction_df.dropna(subset=["ligand", "pdb_id", "residue", "interaction_type", "distance"], inplace=True)

# # Group by Ligand, PDB, Residue, and Chain
# grouped = interaction_df.groupby(["ligand", "pdb_id", "residue", "chain"]).agg({
#     "interaction_type": lambda x: ", ".join(set(x)),  # Unique interaction types
#     "distance": lambda x: ", ".join(set(x))  # Unique distances
# }).reset_index()

# # Format output text
# grouped["output"] = grouped.apply(lambda row: 
#     f"Ligand {row['ligand']} in PDB {row['pdb_id']} interacts with:\n"
#     f" - Residue {row['residue']} (chain {row['chain']}) via {row['interaction_type']} at {row['distance']} Å.", 
#     axis=1
# )

# # Merge interactions per ligand
# final_grouped = grouped.groupby(["ligand", "pdb_id"]).agg({
#     "output": lambda x: "\n".join(x)
# }).reset_index()

# # Generate final structured interaction entries
# aggregated_interaction_data = [
#     {"input": f"What interactions does {row['ligand']} have in PDB {row['pdb_id']}?", "output": row["output"]}
#     for _, row in final_grouped.iterrows()
# ]

# # Filter out old interaction data from the original dataset
# non_interaction_data = df[~df["input"].str.startswith("What interactions does")].to_dict(orient="records")

# # Combine the preserved data with the updated interactions
# final_training_data = non_interaction_data + aggregated_interaction_data

# # Replace any lingering Unicode \u00c5 with proper Å
# final_training_data = replace_unicode_angstrom(final_training_data)

# # Save back to JSON
# with open("training_data.json", "w", encoding="utf-8") as f:
#     json.dump(final_training_data, f, indent=4, ensure_ascii=False)

# print("✅ Training data successfully updated: Interactions aggregated and Angstrom symbol cleaned!")


import json

# Load your data
with open("training_data.json", "r", encoding="utf-8") as f:
    data = json.load(f)

# Clean the data
for entry in data:
    for key in ["input", "output"]:
        if key in entry:
            entry[key] = entry[key].replace(".\n", " ").replace(":\n", ":")

# Save to a new file
with open("training_data_cleaned.json", "w", encoding="utf-8") as f:
    json.dump(data, f, indent=4)

print("Finished cleaning. Saved as training_data_cleaned.json")
