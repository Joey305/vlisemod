import subprocess

# === List your training scripts here (comment out any you want to skip) ===
training_scripts = [
    # "training_1a_generate_data.py",
    # "training_1b_generate_data.py",
    # "training_1c_generate_Contrastivedata.py",
    # "training_1d_generate_Contrastivedata.py",
    # "training_1e_generate_Dialog.py",
    # "training_1eII_generate_Dialog.py",
    # "training_1f_extended.py",
    # "training_1g_Augmented.py",
    # "training_",
    "training_2_Cleaning.py",
    "training_2b_Cleaning.py",
    "training_3Enrich.py",
    "training_4CleanEnrich.py",
    "training_5AdditionalEnrichment.py",
    "training_6CleaningNones.py",
    "training_7Synonyms.py",
    "training_8llama.py",
    # "training_1h_Combination",
    "training_9LlamaFinetuning.py",

]

# === Execute each script one after the other ===
for script in training_scripts:
    print(f"\n🚀 Running {script} ...")
    try:
        subprocess.run(["python", script], check=True)
        print(f"✅ Finished {script}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running {script}: {e}")
        break  # Optional: stop the chain on error
