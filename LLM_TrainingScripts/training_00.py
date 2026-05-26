import os

# Directory with training scripts
directory = "."

# Get all matching scripts in alphanumerical order
training_scripts = sorted([
    f for f in os.listdir(directory)
    if f.startswith("training_") and f.endswith(".py")
])

# Write to catalog file
with open("training_text.txt", "w") as f:
    f.write("# Training Script Catalog\n")
    f.write("# ========================\n\n")

    for script in training_scripts:
        f.write(f"\n### {script}\n")
        f.write("#" * 80 + "\n")
        
        with open(os.path.join(directory, script), "r") as script_file:
            content = script_file.read()
            f.write(content)
            f.write("\n\n" + "-" * 80 + "\n")

print("✅ Full training script catalog written to training_text.txt.")
