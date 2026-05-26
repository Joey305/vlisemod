# import json
# import os

# # --- Helper to convert dialog format to flat "text" format ---
# def flatten_dialog_json(input_file, output_file):
#     with open(input_file) as f:
#         dialogs = json.load(f)

#     flat_data = []

#     for dialog_entry in dialogs:
#         conversation = dialog_entry["dialog"]
#         text = ""
#         for turn in conversation:
#             role = turn["role"]
#             content = turn["content"]
#             if role == "user":
#                 text += f"### User:\n{content}\n\n"
#             elif role == "assistant":
#                 text += f"### Assistant:\n{content}\n\n"
#         flat_data.append({"text": text.strip()})

#     with open(output_file, "w") as f:
#         json.dump(flat_data, f, indent=2)

#     print(f"✅ Flattened {len(flat_data)} dialogs from {input_file} into {output_file}")
#     return flat_data


# # --- Step 1: Convert dialog JSONs ---
# # dialog_short = flatten_dialog_json("training_data_contrastive_dialog.json", "training_data_chat_text.json")
# # dialog_long  = flatten_dialog_json("training_data_contrastive_dialog_long.json", "training_data_chat_text_long.json")

# # --- Step 2: Combine all datasets ---
# combined = []

# # Files already in "input"/"output" format
# standard_files = [
#     "training_data_contrastive.json",
#     "training_data_reasoning_prompts.json",
#     "training_data_multiturn.json",
#     "training_data_chat_text_long.json",
#     "training_data_chat_text.json"
# ]

# for fname in standard_files:
#     if os.path.exists(fname):
#         with open(fname) as f:
#             data = json.load(f)
#             # Convert to "text" format
#             for entry in data:
#                 combined.append({
#                     "text": f"### User:\n{entry['input']}\n\n### Assistant:\n{entry['output']}"
#                 })
#         print(f"✅ Loaded and converted {len(data)} from {fname}")
#     else:
#         print(f"⚠️ Missing file: {fname} — skipping.")

# # Include the newly flattened dialog data
# combined += dialog_short
# combined += dialog_long

# # Save final finetune-ready dataset
# with open("training_data_contrastive_llama.json", "w") as f:
#     json.dump(combined, f, indent=2)

# print(f"\n📦 Final combined dataset: {len(combined)} entries written to llama_finetune_ready.json")



import json
import os

def convert_to_messages_format(entry):
    if "input" in entry and "output" in entry:
        return {
            "messages": [
                {"role": "user", "content": entry["input"]},
                {"role": "assistant", "content": entry["output"]}
            ]
        }
    return None

# === Files to convert ===
input_output_files = [
    "training_data_contrastive.json",
    "training_data_reasoning_prompts.json",
    "training_data_contrastive_multiturn.json"
]

# === Files already in "messages" format ===
chat_files = [
    "training_data_llama.json"
]

# === Files in "text" format that need to be split ===
flattened_text_files = [
    "training_data_chat_text.json",
    "training_data_chat_text_long.json"
]

final_data = []

# Step 1: Convert input/output format
for fname in input_output_files:
    if os.path.exists(fname):
        with open(fname) as f:
            raw = json.load(f)
            converted = [convert_to_messages_format(e) for e in raw if convert_to_messages_format(e)]
            final_data.extend(converted)
            print(f"✅ Converted {len(converted)} entries from {fname}")
    else:
        print(f"⚠️ Missing file: {fname}")

# Step 2: Append raw messages (chat format)
for fname in chat_files:
    if os.path.exists(fname):
        with open(fname) as f:
            chat_data = json.load(f)
            final_data.extend(chat_data)
            print(f"✅ Appended {len(chat_data)} chat entries from {fname}")
    else:
        print(f"⚠️ Missing file: {fname}")

# Step 3: Convert flattened text to messages
for fname in flattened_text_files:
    if os.path.exists(fname):
        with open(fname) as f:
            text_data = json.load(f)
            count = 0
            for item in text_data:
                if "text" in item:
                    parts = item["text"].split("### User:\n")
                    for part in parts:
                        if "### Assistant:\n" in part:
                            user_part, assistant_part = part.split("### Assistant:\n", 1)
                            final_data.append({
                                "messages": [
                                    {"role": "user", "content": user_part.strip()},
                                    {"role": "assistant", "content": assistant_part.strip()}
                                ]
                            })
                            count += 1
            print(f"✅ Converted {count} entries from {fname} (flattened text)")
    else:
        print(f"⚠️ Missing file: {fname}")

# === Save final output ===
# === Save final output safely as JSONL ===
out_path = "training_data_llama_chat_ready.jsonl"
with open(out_path, "w") as f:
    for entry in final_data:
        f.write(json.dumps(entry) + "\n")

print(f"\n📦 Stream-saved {len(final_data)} entries to {out_path}")

print(f"\n📦 Total {len(final_data)} entries written to llama_chat_ready.json")
