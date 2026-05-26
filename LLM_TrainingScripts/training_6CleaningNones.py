import json

# === Load the enriched dataset ===
with open("training_data6.json", "r") as f:
    data = json.load(f)

# === Define filters for bad input patterns ===
def is_bad_entry(entry):
    input_text = entry.get("input", "")
    return (
        input_text.strip() == "Can None be adapted into a PROTAC?" or
        "How similar are ligand None" in input_text
    )

# === Remove bad entries ===
cleaned_data = [entry for entry in data if not is_bad_entry(entry)]
removed_count = len(data) - len(cleaned_data)

# === Overwrite original file with cleaned data ===
with open("training_data7.json", "w") as f:
    json.dump(cleaned_data, f, indent=4)

print(f"✅ Cleaned training_data6.json. Removed {removed_count} bad entries.")
