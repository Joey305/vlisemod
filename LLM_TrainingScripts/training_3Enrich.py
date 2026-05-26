import json
import sqlite3
import re
from copy import deepcopy
import pandas as pd

# === Load synonym and protein mappings ===
conn = sqlite3.connect("viral_data.db")
ligand_synonyms_df = pd.read_sql("SELECT * FROM ligand_synonyms", conn)
virus_proteins_df = pd.read_sql("SELECT * FROM virus_proteins", conn)
conn.close()

# === Build lookup dictionaries ===
ligand_to_synonyms = {}
for _, row in ligand_synonyms_df.iterrows():
    ligand = str(row["ligand"]).strip()
    synonym = str(row["synonym"]).strip()
    if ligand and synonym:  # Only add if both values are non-empty
        ligand_to_synonyms.setdefault(ligand, []).append(synonym)

pdb_to_protein = {}
for _, row in virus_proteins_df.iterrows():
    pdb_id = str(row["pdb_id"]).strip().upper()
    protein = str(row["protein"]).strip()
    if pdb_id and protein:
        pdb_to_protein[pdb_id] = protein

# === Load cleaned training data ===
with open("training_data4.json", "r") as f:
    training_data = json.load(f)

# === Pattern to extract 3-letter ligand and 4-letter PDB code ===
pattern = re.compile(r"What interactions does (\w{3}) have in PDB (\w{4})\?")

# === Build enriched data ===
final_data = []

for entry in training_data:
    input_text = entry["input"]

    # Skip bad entries with 'None'
    if "None" in input_text or "Ligand None" in entry["output"]:
        continue

    match = pattern.match(input_text)
    if match:
        ligand, pdb_id = match.groups()
        synonyms = ligand_to_synonyms.get(ligand, [])
        protein_name = pdb_to_protein.get(pdb_id.upper())

        # ✅ Always include original
        final_data.append(entry)

        # ➕ Synonym variants (ONLY if they exist)
        if synonyms:
            for syn in synonyms:
                syn_entry = deepcopy(entry)
                syn_entry["input"] = f"What interactions does {syn} have in PDB {pdb_id}?"
                syn_entry["output"] = entry["output"].replace(f"Ligand {ligand}", f"Ligand {syn}")
                final_data.append(syn_entry)

        # ➕ Protein name variant
        if protein_name:
            prot_input = f"What interactions does {ligand} have with {protein_name}?"
            prot_output = entry["output"].replace(f"in PDB {pdb_id}", f"with {protein_name}")
            final_data.append({"input": prot_input, "output": prot_output})

            # Combine synonym + protein if synonyms exist
            if synonyms:
                for syn in synonyms:
                    final_data.append({
                        "input": f"What interactions does {syn} have with {protein_name}?",
                        "output": entry["output"]
                            .replace(f"Ligand {ligand}", f"Ligand {syn}")
                            .replace(f"in PDB {pdb_id}", f"with {protein_name}")
                    })
    else:
        final_data.append(entry)

# === Save final JSON ===
with open("training_data5.json", "w") as f:
    json.dump(final_data, f, indent=4)

print("✅ Final training data saved — no 'None' entries, all logic respected.")
