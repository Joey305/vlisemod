import json

# Load the original dataset
with open("training_data8.json", "r") as infile:
    raw_data = json.load(infile)

# Reformat for LLaMA SFT
formatted_data = []
for item in raw_data:
    formatted_data.append({
        "messages": [
            {"role": "user", "content": item["input"]},
            {"role": "assistant", "content": item["output"]}
        ]
    })

# Save to new file
with open("training_data_llama.json", "w") as outfile:
    json.dump(formatted_data, outfile, indent=2)

print("✅ Saved to training_data_llama.json")
