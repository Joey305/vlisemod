import json
import sqlite3
import re
from copy import deepcopy
import pandas as pd

# === Load synonym mappings ===
conn = sqlite3.connect("viral_data.db")
ligand_synonyms_df = pd.read_sql("SELECT * FROM ligand_synonyms", conn)
conn.close()

# === Build synonym dictionary ===
ligand_to_synonyms = {}
for _, row in ligand_synonyms_df.iterrows():
    ligand = str(row["ligand"]).strip()
    synonym = str(row["synonym"]).strip()
    if ligand and synonym:
        ligand_to_synonyms.setdefault(ligand, []).append(synonym)

# === Load the existing training dataset ===
with open("training_data5.json", "r") as f:
    original_data = json.load(f)

# === Regex patterns ===
similarity_pattern = re.compile(r"How similar are ligand (\w{3}) and ligand (\w{3})\?")
protac_pattern = re.compile(r"Can (\w{3}) be adapted into a PROTAC\?")

# === Begin enrichment ===
enriched_data = []
for entry in original_data:
    enriched_data.append(entry)  # Always keep original
    input_text = entry["input"]

    # 🔬 Similarity enrichment
    match_sim = similarity_pattern.match(input_text)
    if match_sim:
        lig1, lig2 = match_sim.groups()
        syns1 = ligand_to_synonyms.get(lig1, [])
        syns2 = ligand_to_synonyms.get(lig2, [])

        for s1 in syns1:
            new_entry = deepcopy(entry)
            new_entry["input"] = f"How similar are ligand {s1} and ligand {lig2}?"
            new_entry["output"] = entry["output"].replace(f"Ligand {lig1}", f"Ligand {s1}")
            enriched_data.append(new_entry)

        for s2 in syns2:
            new_entry = deepcopy(entry)
            new_entry["input"] = f"How similar are ligand {lig1} and ligand {s2}?"
            new_entry["output"] = entry["output"].replace(f"Ligand {lig2}", f"Ligand {s2}")
            enriched_data.append(new_entry)

        for s1 in syns1:
            for s2 in syns2:
                new_entry = deepcopy(entry)
                new_entry["input"] = f"How similar are ligand {s1} and ligand {s2}?"
                new_entry["output"] = entry["output"].replace(f"Ligand {lig1}", f"Ligand {s1}").replace(f"Ligand {lig2}", f"Ligand {s2}")
                enriched_data.append(new_entry)

    # 🔧 PROTAC enrichment
    match_protac = protac_pattern.match(input_text)
    if match_protac:
        lig = match_protac.group(1)
        synonyms = ligand_to_synonyms.get(lig, [])
        for syn in synonyms:
            new_entry = deepcopy(entry)
            new_entry["input"] = f"Can {syn} be adapted into a PROTAC?"
            new_entry["output"] = entry["output"].replace(f"PROTAC Adaptability Score for {lig}", f"PROTAC Adaptability Score for {syn}")
            enriched_data.append(new_entry)

# === Save enriched version ===
with open("training_data6.json", "w") as f:
    json.dump(enriched_data, f, indent=4)

print(f"✅ Enriched dataset saved to training_data6.json with {len(enriched_data) - len(original_data)} new entries.")
