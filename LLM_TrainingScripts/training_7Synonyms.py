import sqlite3
import json
import re
from copy import deepcopy

# === Load training data ===
with open("training_data7.json", "r") as f:
    data = json.load(f)

# === Connect to DB and load mappings ===
conn = sqlite3.connect("viral_data.db")
cursor = conn.cursor()

# PDB → Protein mapping
cursor.execute("SELECT pdb_id, protein FROM Virus_Proteins")
pdb_to_protein = {row[0].upper(): row[1] for row in cursor.fetchall()}

# Ligand → Synonyms mapping (safe check for None)
cursor.execute("SELECT ligand, synonym FROM ligand_synonyms")
ligand_to_synonyms = {}
for ligand, synonym in cursor.fetchall():
    if ligand and synonym:
        ligand = ligand.strip().upper()
        synonym = synonym.strip()
        ligand_to_synonyms.setdefault(ligand, []).append(synonym)

conn.close()

# Regex to match "What interactions does XXX (also known as ...) have in PDB XXXX?"
pattern = re.compile(r"What interactions does (\w{3}) \(also known as .+?\) have in PDB (\w{4})\?", re.IGNORECASE)

final_data = []
new_entries = 0

for entry in data:
    input_text = entry["input"]
    output_text = entry["output"]
    final_data.append(entry)  # Always keep the original

    match = pattern.match(input_text)
    if match:
        ligand, pdb_id = match.groups()
        ligand = ligand.upper()
        pdb_id = pdb_id.upper()

        protein = pdb_to_protein.get(pdb_id)
        synonyms = ligand_to_synonyms.get(ligand, [])

        if protein:
            # ➕ Ligand with protein
            final_data.append({
                "input": f"What interactions does {ligand} have with {protein}?",
                "output": output_text.replace(f"in PDB {pdb_id}", f"with {protein}")
            })
            new_entries += 1

            # ➕ Each synonym with protein
            for syn in synonyms:
                final_data.append({
                    "input": f"What interactions does {syn} have with {protein}?",
                    "output": output_text
                        .replace(f"Ligand {ligand}", f"Ligand {syn}")
                        .replace(f"in PDB {pdb_id}", f"with {protein}")
                })
                new_entries += 1

# === Save output ===
with open("training_data8.json", "w") as f:
    json.dump(final_data, f, indent=2)

print(f"✅ Enriched with {new_entries} new entries using protein + synonym logic.")
print(f"📦 Saved to training_data8.json with total {len(final_data)} entries.")
