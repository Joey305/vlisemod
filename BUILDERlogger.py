import os
import re
import csv
from rdkit import Chem
from datetime import datetime

LOG_DIR = "copy_logs"
CSV_PATH = "static/data/PROTAC_log.csv"


# -----------------------------
# Helper: Convert MOL block to SMART wildcard SMILES
# -----------------------------
def molblock_to_smiles(molblock):
    if not molblock.strip():
        return ""

    mol = Chem.MolFromMolBlock(molblock, sanitize=False, removeHs=False)
    if mol is None:
        print("⚠️ RDKit could not load MOL block.")
        return ""

    for atom in mol.GetAtoms():
        label = atom.GetProp("_MolFileAtomLabel") if atom.HasProp("_MolFileAtomLabel") else ""
        if label == "R1":
            atom.SetAtomicNum(0)
            atom.SetAtomMapNum(1)
        elif label == "R2":
            atom.SetAtomicNum(0)
            atom.SetAtomMapNum(2)

    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except:
        return ""


# -----------------------------
# Helper: extract MOL sections
# -----------------------------
def extract_section(text, label):
    pattern = rf"=== {label} ===\n(.*?)(?====|\Z)"
    match = re.search(pattern, text, flags=re.S)
    return match.group(1).strip() if match else ""


# -----------------------------
# Process a single .txt log into dict
# -----------------------------
def process_log_file(path):
    with open(path, "r") as f:
        txt = f.read()

    ip_match = re.search(r"IP:\s*(.*)", txt)
    ip = ip_match.group(1).strip() if ip_match else "UNKNOWN"

    warhead_mol = extract_section(txt, "WARHEAD")
    linker_mol = extract_section(txt, "LINKER")
    ligase_mol = extract_section(txt, "LIGASE")

    protac_section = extract_section(txt, "PROTAC \\(MOL \\+ SMILES\\)")
    protac_mol = "\n".join(protac_section.split("\n")[:-1])
    protac_smiles = protac_section.split("\n")[-1].strip()

    warhead_smiles = molblock_to_smiles(warhead_mol)
    linker_smiles = molblock_to_smiles(linker_mol)
    ligase_smiles = molblock_to_smiles(ligase_mol)
    protac_smiles_rdkit = molblock_to_smiles(protac_mol)

    return {
        "date": os.path.basename(path).replace("protac_", "").replace(".txt", ""),
        "ip": ip,
        "warhead": warhead_smiles,
        "linker": linker_smiles,
        "ligase": ligase_smiles,
        "protac": protac_smiles_rdkit or protac_smiles
    }


# -----------------------------
# Read latest timestamp already in CSV
# -----------------------------
def get_latest_csv_timestamp():
    if not os.path.exists(CSV_PATH):
        return None

    with open(CSV_PATH, "r") as f:
        rows = list(csv.reader(f))

    if len(rows) <= 1:
        return None

    # Last row has latest date because CSV is append-only
    last_date = rows[-1][0]

    try:
        return datetime.strptime(last_date, "%Y-%m-%d_%H-%M-%S")
    except:
        return None


# -----------------------------
# Append a row to the CSV
# -----------------------------
def append_to_csv(row):
    with open(CSV_PATH, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            row["date"], row["ip"], row["protac"],
            row["warhead"], row["linker"], row["ligase"]
        ])


# -----------------------------
# Ensures header exists
# -----------------------------
def ensure_csv_header(path):
    if not os.path.exists(path):
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Date", "IP", "PROTAC", "WARHEAD", "LINKER", "LIGASE"])


# -----------------------------
# Main: only append new logs
# -----------------------------
def main():
    ensure_csv_header(CSV_PATH)
    last_ts = get_latest_csv_timestamp()

    print(f"📅 Latest timestamp in CSV: {last_ts}\n")

    files = sorted([f for f in os.listdir(LOG_DIR) if f.endswith(".txt")])
    print(f"📂 Found {len(files)} .txt log files.\n")

    for fname in files:
        date_str = fname.replace("protac_", "").replace(".txt", "")
        file_ts = datetime.strptime(date_str, "%Y-%m-%d_%H-%M-%S")

        if last_ts and file_ts <= last_ts:
            print(f"⏭️ Skipping {fname} (already logged)")
            continue

        print(f"➕ Adding NEW log: {fname}")
        row = process_log_file(os.path.join(LOG_DIR, fname))
        append_to_csv(row)

    print("\n🎉 DONE — Only new logs appended!")


if __name__ == "__main__":
    main()
