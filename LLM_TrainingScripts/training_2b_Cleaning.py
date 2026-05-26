import json

# Load your data
with open("training_data3.json", "r", encoding="utf-8") as f:
    data = json.load(f)

# Clean the data
for entry in data:
    for key in ["input", "output"]:
        if key in entry:
            entry[key] = entry[key].replace(".\n", " ").replace(":\n", ":")

# Save to a new file
with open("training_data4.json", "w", encoding="utf-8") as f:
    json.dump(data, f, indent=4)

print("Finished cleaning. Saved as training_data4.json")


