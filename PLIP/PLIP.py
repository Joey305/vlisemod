# import os
# from plip.structure.preparation import PDBComplex

# # Get user input for PDB file, ligand, chain, and residue number
# pdb_file = input("Enter the PDB file name (e.g., 3eky.pdb): ").strip()
# ligand_name = input("Enter the ligand name (e.g., DR7): ").strip().upper()
# chain_id = input("Enter the ligand chain (e.g., A): ").strip().upper()
# residue_number = input("Enter the ligand residue number (e.g., 100): ").strip()

# # Validate input
# if not os.path.exists(pdb_file):
#     print(f"Error: PDB file '{pdb_file}' not found.")
#     exit(1)

# # Standardized naming for output
# session_name = f"{pdb_file.split('.')[0]}_{ligand_name}_{residue_number}{chain_id}"
# output_pml = f"{session_name}.pml"

# # Load PDB into PLIP
# print("\n[INFO] Loading PDB into PLIP...")
# my_mol = PDBComplex()
# my_mol.load_pdb(pdb_file) 
# my_mol.analyze()

# # Construct ligand ID (HetID:Chain:ResidueNumber)
# ligand_id = f"{ligand_name}:{chain_id}:{residue_number}"

# # Extract interaction data
# if ligand_id in my_mol.interaction_sets:
#     my_interactions = my_mol.interaction_sets[ligand_id]

#     # Print interactions of interest
#     print("\n[INFO] Ligand-Protein Interactions Found:")
#     print("Residues involved in pi-stacking:", [pistack.resnr for pistack in my_interactions.pistacking])
#     print("Residues involved in hydrophobic interactions:", [hydrophobic.resnr for hydrophobic in my_interactions.hydrophobic_contacts])

#     # Generate PyMOL script
#     print(f"\n[INFO] Generating PyMOL script: {output_pml}")

#     with open(output_pml, "w") as pymol_script:
#         pymol_script.write(f"""
# delete all
# load {pdb_file}
# select ligand, resn {ligand_name} and chain {chain_id} and resi {residue_number}
# select protein, polymer.protein
# show cartoon, protein
# show sticks, ligand
# color yellow, ligand
# show lines, ligand
# """)
    
#     print(f"\n✅ Analysis complete! Open in PyMOL with:\n  pymol -r {output_pml}")

# else:
#     print("\n❌ Error: Ligand binding site not found in the PDB file.")
#     print("Check if the residue number, chain, or ligand ID is correct.")






from plip.structure.preparation import PDBComplex
from plip.visualization.pymol import PyMOLVisualizer

# Load the PDB file into PLIP
my_mol = PDBComplex()
my_mol.load_pdb('3eky.pdb')  # Ensure the PDB file is in the working directory

# Analyze interactions
my_mol.analyze()

# Unique binding site identifier (HetID:Chain:Position)
my_bsid = 'DR7:A:100'

# Ensure the ligand binding site exists before proceeding
if my_bsid in my_mol.interaction_sets:
    my_interactions = my_mol.interaction_sets[my_bsid]  # Extract interaction data

    # ✅ PRINT ALL INTERACTING RESIDUES
    interacting_residues = set()
    
    # Collect all interacting residue numbers from all interactions
    for interaction in my_interactions.all_itypes:
        for inter in my_interactions.interactions[interaction]:  # Access each interaction list
            interacting_residues.add(inter.resnr)  # Store residue numbers

    print("All interacting residues:", sorted(interacting_residues))

    # ✅ FIXED: Initialize PyMOL Visualizer with only `PDBComplex`
    visualizer = PyMOLVisualizer(my_mol)

    # Generate PyMOL session for the DR7 ligand
    visualizer.visualize_site(my_interactions)  # Creates PyMOL visualization
    visualizer.save_session("3eky_DR7.pse")  # Saves as PyMOL session file

    print("✅ PyMOL session file generated: 3eky_DR7.pse")
else:
    print(f"⚠️ No interactions found for {my_bsid}")
