# import pandas as pd
# import json
# import sqlite3
# from rdkit import Chem
# from rdkit.Chem import DataStructs, Descriptors, QED
# from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

# # Connect to the SQLite database
# conn = sqlite3.connect("viral_data.db")

# # Load tables dynamically
# tables = {
#     "ligand_synonyms": "ligand_synonyms",
#     "functional_grouped": "Functional_GROUPED",
#     "arpeggio_contacts": "Arpeggio_Contacts_Data",
#     "virus_proteins": "Virus_Proteins",
#     "ligand_atoms": "ligand_atoms",
#     "rupley_sasa": "RUPLEY_SASA_DATA",
#     "ligand_water_distances": "ligand_water_distances",
#     "receptor_binding_pocket": "receptor_binding_pocket",
#     "ligand_atoms_smiles": "Ligand_Atoms_Smiles",
#     "binding_affinity_data": "binding_affinity"
# }

# # Dictionary to hold loaded tables
# dataframes = {}

# for key, query in tables.items():
#     try:
#         dataframes[key] = pd.read_sql_query(f"SELECT * FROM {query}", conn)
#     except Exception as e:
#         print(f"⚠️ Warning: Could not load {query}. Skipping...")
#         dataframes[key] = None  # Assign None if the table is missing

# # Close database connection
# conn.close()

# # Unpack dataframes
# ligand_synonyms = dataframes["ligand_synonyms"]
# functional_grouped = dataframes["functional_grouped"]
# arpeggio_contacts = dataframes["arpeggio_contacts"]
# virus_proteins = dataframes["virus_proteins"]
# ligand_atoms = dataframes["ligand_atoms"]
# rupley_sasa = dataframes["rupley_sasa"]
# ligand_water_distances = dataframes["ligand_water_distances"]
# receptor_binding_pocket = dataframes["receptor_binding_pocket"]
# ligand_atoms_smiles = dataframes["ligand_atoms_smiles"]
# binding_affinity_data = dataframes["binding_affinity_data"]

# # Synonym lookup dictionary
# synonym_lookup = {}
# if ligand_synonyms is not None:
#     for _, row in ligand_synonyms.iterrows():
#         if pd.isna(row.get("ligand")) or pd.isna(row.get("synonym")):
#             continue  # Skip missing values
#         synonym_lookup.setdefault(row["ligand"], []).append(row["synonym"])

# def get_ligand_name(ligand):
#     """Get the ligand name with synonyms."""
#     synonyms = synonym_lookup.get(ligand, [])
#     return f"{ligand} (also known as {', '.join(synonyms)})" if synonyms else ligand

# # Training data storage
# training_data = []

# # Generate training data
# datasets = [
#     ("ligand_synonyms", "ligand", "synonym", "What is another name for {ligand}?", "{ligand} is also known as {synonym}."),
#     ("functional_grouped", "ligand", "functional_groups", "What functional groups are in {ligand_name} in PDB {pdb_id}?",
#      "Ligand {ligand_name} in PDB {pdb_id} contains: {functional_groups}."),
#     ("arpeggio_contacts", "ligand", "residue", "What interactions does {ligand_name} have in PDB {pdb_id}?",
#      "Ligand {ligand_name} in PDB {pdb_id} interacts with residue {residue} {residue_number} (chain {residue_chain}) via a {Contact} at {Distance} Å."),
#     ("rupley_sasa", "ligand", "SASA_Area", "What is the SASA value for atom {atom_id} in {ligand_name} in PDB {pdb_id}?",
#      "The SASA value for atom {atom_id} in {ligand_name} in PDB {pdb_id} is {SASA_Area} Å²."),
#     ("ligand_water_distances", "ligand", "distance", "What is the distance between {ligand_name} atom {atom_id} and water {water_sequence_id} in PDB {pdb_id}?",
#      "The distance between {ligand_name} atom {atom_id} and water {water_sequence_id} in PDB {pdb_id} is {distance} Å."),
#     ("receptor_binding_pocket", "residue", "pdb_id", "What residues are in the binding pocket of PDB {pdb_id}?",
#      "The binding pocket of PDB {pdb_id} includes residue {residue} {residue_number} (chain {residue_chain})."),
# ]

# for dataset, key, col, input_templ, output_templ in datasets:
#     df = dataframes[dataset]
#     if df is not None:
#         for _, row in df.iterrows():
#             try:
#                 if pd.isna(row.get(key)) or pd.isna(row.get(col)):
#                     continue  # Skip missing values
#                 ligand_name = get_ligand_name(row[key]) if key == "ligand" else row[key]
#                 training_data.append({"input": input_templ.format(**row, ligand_name=ligand_name),
#                                       "output": output_templ.format(**row, ligand_name=ligand_name)})
#             except Exception as e:
#                 print(f"⚠️ Error processing {dataset}: {e}")

# # RDKit Ligand Similarity Calculation
# if ligand_atoms_smiles is not None:
#     valid_smiles_df = ligand_atoms_smiles.dropna(subset=["smiles"])
#     valid_smiles_df = valid_smiles_df[valid_smiles_df["smiles"].apply(lambda x: Chem.MolFromSmiles(x) is not None)]
#     valid_ligands = valid_smiles_df["ligand"].unique()

#     def get_valid_molecule(smiles):
#         try:
#             mol = Chem.MolFromSmiles(smiles)
#             Chem.SanitizeMol(mol)
#             return mol
#         except:
#             return None

#     for i, ligand1 in enumerate(valid_ligands):
#         for ligand2 in valid_ligands[i+1:]:
#             try:
#                 mol1, mol2 = get_valid_molecule(valid_smiles_df[valid_smiles_df["ligand"] == ligand1]["smiles"].values[0]), \
#                              get_valid_molecule(valid_smiles_df[valid_smiles_df["ligand"] == ligand2]["smiles"].values[0])

#                 if mol1 and mol2:
#                     generator = GetMorganGenerator(radius=2)
#                     sim_score = DataStructs.TanimotoSimilarity(generator.GetFingerprint(mol1), generator.GetFingerprint(mol2))

#                     training_data.append({
#                         "input": f"How similar are ligand {ligand1} and ligand {ligand2}?",
#                         "output": f"Ligand {ligand1} and ligand {ligand2} have a Tanimoto similarity score of {sim_score:.2f}."
#                     })
#             except Exception as e:
#                 print(f"⚠️ Error processing similarity for {ligand1}, {ligand2}: {e}")

# # Binding affinity data processing
# if binding_affinity_data is not None:
#     for _, row in binding_affinity_data.iterrows():
#         try:
#             if pd.isna(row.get("ligand")) or pd.isna(row.get("protein")) or pd.isna(row.get("binding_affinity")):
#                 continue  # Skip missing values

#             ligand_name = get_ligand_name(row["ligand"])
#             training_data.append({
#                 "input": f"What is the binding affinity of ligand {ligand_name} with protein {row['protein']}?",
#                 "output": f"Ligand {ligand_name} binds to protein {row['protein']} with an affinity of {row['binding_affinity']} kcal/mol."
#             })
#         except Exception as e:
#             print(f"⚠️ Error processing binding affinity: {e}")

# # Save training data to JSON
# with open("training_data.json", "w") as f:
#     json.dump(training_data, f, indent=4)

# print("✅ Training data generated and saved to training_data.json!")





import pandas as pd
import json
import sqlite3
import re
from rdkit import Chem
from rdkit.Chem import Descriptors

# Connect to the SQLite database
conn = sqlite3.connect("viral_data.db")

# Load necessary tables
try:
    functional_group_atoms = pd.read_sql_query("SELECT * FROM Functional_Group_Atoms", conn)
    rupley_sasa = pd.read_sql_query("SELECT * FROM RUPLEY_SASA_DATA", conn)
    ligand_smiles = pd.read_sql_query("SELECT * FROM Ligand_Atoms_Smiles", conn)
except Exception as e:
    print(f"⚠️ Warning: Could not load some tables: {e}")
    conn.close()
    exit()

# Close the database connection
conn.close()

# Load existing training data
try:
    with open("training_data.json", "r") as f:
        training_data = json.load(f)
except FileNotFoundError:
    training_data = []

# List of functional groups favorable for PROTAC linkers
PROTAC_LINKER_GROUPS = ["Hydroxyl/Alcohol", "Amine", "Carboxylic-Acid", "Thiol"]

# Function to clean functional group names by removing Roman numerals
def clean_functional_group_name(name):
    return re.sub(r"\s*\(.*?\)", "", name).strip()

# Update the functional group column to remove Roman numerals
functional_group_atoms["functional_group"] = functional_group_atoms["functional_group"].apply(clean_functional_group_name)

# Function to calculate molecular weight from SMILES
def calculate_molecular_weight(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Descriptors.MolWt(mol)
    except Exception as e:
        print(f"⚠️ Error calculating molecular weight for SMILES {smiles}: {e}")
    return None

# Add molecular weight data to the Ligand_Atoms_Smiles table
ligand_smiles["molecular_weight"] = ligand_smiles["smiles"].apply(calculate_molecular_weight)

# Function to assess PROTAC potential
def assess_protac_potential(ligand):
    """Determine if a ligand is suitable for PROTAC development based on solvent-exposed functional groups and molecular weight."""
    # Get functional group data for the ligand
    functional_groups = functional_group_atoms[functional_group_atoms["ligand"] == ligand]

    if functional_groups.empty:
        return 0, "No functional group data available."

    # Get SASA data for the ligand
    sasa_data = rupley_sasa[rupley_sasa["ligand"] == ligand]

    if sasa_data.empty:
        return 0, "No solvent accessibility data available."

    # Get molecular weight for the ligand
    smiles_row = ligand_smiles[ligand_smiles["ligand"] == ligand]
    molecular_weight = smiles_row["molecular_weight"].values[0] if not smiles_row.empty else None

    if molecular_weight is None:
        return 0, "No molecular weight data available."

    # Find solvent-exposed functional group atoms
    exposed_groups = []
    for _, row in functional_groups.iterrows():
        atom_id = row["atom_id"]
        group = row["functional_group"]

        if group in PROTAC_LINKER_GROUPS:
            sasa_entry = sasa_data[sasa_data["atom_id"] == atom_id]
            
            # Check if the atom is solvent-exposed
            if not sasa_entry.empty:
                exposed_groups.append(group)

    # Assess PROTAC potential based on the number of solvent-exposed functional groups and molecular weight
    if molecular_weight > 500:
        weight_comment = f"However, its molecular weight ({molecular_weight:.2f} Da) is higher than the optimal range for PROTAC development."
    else:
        weight_comment = f"Its molecular weight ({molecular_weight:.2f} Da) is within the optimal range for PROTAC development."

    if len(exposed_groups) > 1:
        score = 3
        reasoning = f"This compound has multiple exposed linker-friendly functional groups: {', '.join(set(exposed_groups))}. {weight_comment}"
    elif len(exposed_groups) == 1:
        score = 2
        reasoning = f"This compound has an exposed {exposed_groups[0]} functional group, making it a viable PROTAC candidate. {weight_comment}"
    else:
        score = 1
        reasoning = f"This compound has no solvent-exposed linker-friendly groups and may not be ideal for PROTAC development. {weight_comment}"

    return score, reasoning

# Generate PROTAC-related training data
for ligand in functional_group_atoms["ligand"].unique():
    try:
        protac_score, reasoning = assess_protac_potential(ligand)

        input_text = f"Can {ligand} be adapted into a PROTAC?"
        output_text = f"PROTAC Adaptability Score for {ligand}: {protac_score}/3. {reasoning}"

        training_data.append({"input": input_text, "output": output_text})

    except Exception as e:
        print(f"⚠️ Error processing {ligand}: {e}")

# Save updated training data to JSON
with open("training_data.json", "w") as f:
    json.dump(training_data, f, indent=4)

print("✅ PROTAC Adaptability Assessment complete. Data appended to training_data.json!")
