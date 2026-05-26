import json

# Path to your training data file
json_path = "training_data5.json"

# Load the JSON
with open(json_path, "r") as f:
    data = json.load(f)

# Filter out entries with "What interactions does None" in the input
cleaned_data = [entry for entry in data if not entry["input"].startswith("What interactions does None")]

# Save it back to the same file
with open(json_path, "w") as f:
    json.dump(cleaned_data, f, indent=4)

print(f"✅ Cleaned {len(data) - len(cleaned_data)} 'None' entries. Saved to {json_path}")
