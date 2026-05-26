# import pandas as pd
# import json
# import sqlite3
# import re
# from rdkit import Chem
# from rdkit.Chem import Descriptors

# # Connect to the SQLite database
# conn = sqlite3.connect("viral_data.db")

# # Load necessary tables
# try:
#     functional_group_atoms = pd.read_sql_query("SELECT * FROM Functional_Group_Atoms", conn)
#     rupley_sasa = pd.read_sql_query("SELECT * FROM RUPLEY_SASA_DATA", conn)
#     ligand_smiles = pd.read_sql_query("SELECT * FROM Ligand_Atoms_Smiles", conn)
# except Exception as e:
#     print(f"⚠️ Warning: Could not load some tables: {e}")
#     conn.close()
#     exit()

# # Close the database connection
# conn.close()

# # Load existing training data
# try:
#     with open("training_data1.json", "r") as f:
#         training_data = json.load(f)
# except FileNotFoundError:
#     training_data = []

# # List of functional groups favorable for PROTAC linkers
# PROTAC_LINKER_GROUPS = ["Hydroxyl/Alcohol", "Amine", "Carboxylic-Acid", "Thiol"]

# # Function to clean functional group names by removing Roman numerals
# def clean_functional_group_name(name):
#     return re.sub(r"\s*\(.*?\)", "", name).strip()

# # Update the functional group column to remove Roman numerals
# functional_group_atoms["functional_group"] = functional_group_atoms["functional_group"].apply(clean_functional_group_name)

# # Function to calculate molecular weight from SMILES
# def calculate_molecular_weight(smiles):
#     try:
#         mol = Chem.MolFromSmiles(smiles)
#         if mol:
#             return Descriptors.MolWt(mol)
#     except Exception as e:
#         print(f"⚠️ Error calculating molecular weight for SMILES {smiles}: {e}")
#     return None

# # Add molecular weight data to the Ligand_Atoms_Smiles table
# ligand_smiles["molecular_weight"] = ligand_smiles["smiles"].apply(calculate_molecular_weight)

# # Function to assess PROTAC potential
# def assess_protac_potential(ligand):
#     """Determine if a ligand is suitable for PROTAC development based on solvent-exposed functional groups and molecular weight."""
#     # Get functional group data for the ligand
#     functional_groups = functional_group_atoms[functional_group_atoms["ligand"] == ligand]

#     if functional_groups.empty:
#         return 0, "No functional group data available."

#     # Get SASA data for the ligand
#     sasa_data = rupley_sasa[rupley_sasa["ligand"] == ligand]

#     if sasa_data.empty:
#         return 0, "No solvent accessibility data available."

#     # Get molecular weight for the ligand
#     smiles_row = ligand_smiles[ligand_smiles["ligand"] == ligand]
#     molecular_weight = smiles_row["molecular_weight"].values[0] if not smiles_row.empty else None

#     if molecular_weight is None:
#         return 0, "No molecular weight data available."

#     # Find solvent-exposed functional group atoms
#     exposed_groups = []
#     for _, row in functional_groups.iterrows():
#         atom_id = row["atom_id"]
#         group = row["functional_group"]

#         if group in PROTAC_LINKER_GROUPS:
#             sasa_entry = sasa_data[sasa_data["atom_id"] == atom_id]
            
#             # Check if the atom is solvent-exposed
#             if not sasa_entry.empty:
#                 exposed_groups.append(group)

#     # Assess PROTAC potential based on the number of solvent-exposed functional groups and molecular weight
#     if molecular_weight > 500:
#         weight_comment = f"However, its molecular weight ({molecular_weight:.2f} Da) is higher than the optimal range for PROTAC development."
#     else:
#         weight_comment = f"Its molecular weight ({molecular_weight:.2f} Da) is within the optimal range for PROTAC development."

#     if len(exposed_groups) > 1:
#         score = 3
#         reasoning = f"This compound has multiple exposed linker-friendly functional groups: {', '.join(set(exposed_groups))}. {weight_comment}"
#     elif len(exposed_groups) == 1:
#         score = 2
#         reasoning = f"This compound has an exposed {exposed_groups[0]} functional group, making it a viable PROTAC candidate. {weight_comment}"
#     else:
#         score = 1
#         reasoning = f"This compound has no solvent-exposed linker-friendly groups and may not be ideal for PROTAC development. {weight_comment}"

#     return score, reasoning

# # Generate PROTAC-related training data
# for ligand in functional_group_atoms["ligand"].unique():
#     try:
#         protac_score, reasoning = assess_protac_potential(ligand)

#         input_text = f"Can {ligand} be adapted into a PROTAC?"
#         output_text = f"PROTAC Adaptability Score for {ligand}: {protac_score}/3. {reasoning}"

#         training_data.append({"input": input_text, "output": output_text})

#     except Exception as e:
#         print(f"⚠️ Error processing {ligand}: {e}")

# # Save updated training data to JSON
# with open("training_data2.json", "w") as f:
#     json.dump(training_data, f, indent=4)

# print("✅ PROTAC Adaptability Assessment complete. Data appended to training_data.json!")




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
    arpeggio_contacts = pd.read_sql_query("SELECT * FROM Arpeggio_Contacts_Data", conn)
except Exception as e:
    print(f"⚠️ Warning: Could not load some tables: {e}")
    conn.close()
    exit()

# Close the database connection
conn.close()

# Load existing training data
try:
    with open("training_data1.json", "r") as f:
        training_data = json.load(f)
except FileNotFoundError:
    training_data = []

# List of functional groups favorable for PROTAC linkers
PROTAC_LINKER_GROUPS = ["Hydroxyl/Alcohol", "Amine", "Carboxylic-Acid", "Thiol"]

# Define what counts as a 'high-quality' interaction
# Define what counts as a 'high-quality' interaction
HIGH_QUALITY_CONTACTS = {
    "polar",
    "weak_polar",
    "hbond",
    "ionic",
    "weak_hbond",
    "vdw_clash"
}

successful_ligands = []
failed_ligands = []

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

# Compute binding interaction score
# Compute binding interaction score
def compute_interaction_score(ligand):
    if arpeggio_contacts.empty:
        return 0.0, "No interaction data available.", 0.0, 0

    # Filter for the ligand and high-quality contacts
    filtered = arpeggio_contacts[
        (arpeggio_contacts["ligand"] == ligand) &
        (arpeggio_contacts["Contact"].isin(HIGH_QUALITY_CONTACTS))
    ]

    if filtered.empty:
        return 0.0, "No high-quality contacts found.", 0.0, 0

    # Group by (pdb_id, ligand_id, chain) to avoid counting duplicates
    grouped = filtered.groupby(["pdb_id", "ligand_id", "chain"])
    contact_counts = grouped.size()

    # Take the average number of contacts across unique groups
    avg_contacts = contact_counts.mean()
    n_instances = len(contact_counts)

    # Now scale into score buckets
    if avg_contacts >= 5:
        score = 1.0
        msg = "Strong binder"
    elif avg_contacts >= 2:
        score = 0.5
        msg = "Moderate binder"
    else:
        score = 0.0
        msg = "Weak binder"

    reasoning = f"{msg}: {avg_contacts:.1f} average high-quality contacts across {n_instances} structure(s)."
    return score, reasoning, avg_contacts, n_instances



# Function to assess PROTAC potential
# Function to assess PROTAC potential
def assess_protac_potential(ligand):
    functional_groups = functional_group_atoms[functional_group_atoms["ligand"] == ligand]
    if functional_groups.empty:
        return 0, "No functional group data available.", True

    sasa_data = rupley_sasa[rupley_sasa["ligand"] == ligand]
    if sasa_data.empty:
        return 0, "No solvent accessibility data available.", True

    smiles_row = ligand_smiles[ligand_smiles["ligand"] == ligand]
    molecular_weight = smiles_row["molecular_weight"].values[0] if not smiles_row.empty else None
    if molecular_weight is None:
        return 0, "No molecular weight data available.", True

    exposed_groups = []
    for _, row in functional_groups.iterrows():
        atom_id = row["atom_id"]
        group = row["functional_group"]
        if group in PROTAC_LINKER_GROUPS:
            sasa_entry = sasa_data[sasa_data["atom_id"] == atom_id]
            if not sasa_entry.empty:
                exposed_groups.append(group)

    if molecular_weight > 500:
        weight_comment = f"However, its molecular weight ({molecular_weight:.2f} Da) is higher than the optimal range for PROTAC development."
    else:
        weight_comment = f"Its molecular weight ({molecular_weight:.2f} Da) is within the optimal range for PROTAC development."

    if len(exposed_groups) > 1:
        base_score = 3
        reasoning = f"Multiple exposed linker-friendly functional groups: {', '.join(set(exposed_groups))}. {weight_comment}"
    elif len(exposed_groups) == 1:
        base_score = 2
        reasoning = f"One exposed {exposed_groups[0]} group, a viable PROTAC candidate. {weight_comment}"
    else:
        base_score = 1
        reasoning = f"No solvent-exposed linker-friendly groups. May not be ideal. {weight_comment}"

    interaction_score, interaction_reasoning, avg_contacts, n_structures = compute_interaction_score(ligand)
    interaction_density = avg_contacts / molecular_weight if molecular_weight else 0

    combined_reasoning = (
        f"{reasoning}\n"
        f"Binding Interaction Score: {interaction_score}/1. {interaction_reasoning}\n"
        f"Interaction Density: {interaction_density:.4f} high-quality contacts per Dalton."
    )

    total_score = base_score + interaction_score
    return total_score, combined_reasoning, False, avg_contacts, interaction_density



# Generate PROTAC-related training data
for ligand in functional_group_atoms["ligand"].unique():
    try:
        total_score, reasoning, failed, avg_contacts, interaction_density = assess_protac_potential(ligand)


        if failed:
            # Grab PDB ID for context if available
            try:
                example_row = functional_group_atoms[functional_group_atoms["ligand"] == ligand].iloc[0]
                pdb_id = example_row["pdb_id"]
            except:
                pdb_id = "N/A"

            failed_ligands.append({
                "ligand": ligand,
                "pdb_id": pdb_id,
                "failure_reason": reasoning
            })
            continue  # Skip training data append

        # Only executed if not failed
        example_row = functional_group_atoms[functional_group_atoms["ligand"] == ligand].iloc[0]
        pdb_id = example_row["pdb_id"]

        successful_ligands.append({
            "ligand": ligand,
            "pdb_id": pdb_id,
            "score": total_score,
            "avg_contacts": avg_contacts,
            "interaction_density": interaction_density,
            "reasoning": reasoning
        })


        training_data.append({
            "input": f"Can {ligand} be adapted into a PROTAC?",
            "output": f"PROTAC Adaptability Score for {ligand}: {total_score}/4. {reasoning}"
        })

    except Exception as e:
        failed_ligands.append({
            "ligand": ligand,
            "pdb_id": "N/A",
            "failure_reason": f"Unhandled exception: {e}"
        })



# Save updated training data to JSON
with open("training_data2.json", "w") as f:
    json.dump(training_data, f, indent=4)
# Export successful ligand-PDB-score mapping
# Save failed ligands and reasons
fail_df = pd.DataFrame(failed_ligands)
fail_df.to_csv("failed_protac_candidates.csv", index=False)
print("🧪 Logged failed candidates to failed_protac_candidates.csv")
