#!/usr/bin/env python3
import os
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

# --- CONFIG ---
LIGASE_DIR = "Ligases"
OUT_DIR = os.path.join(LIGASE_DIR, "Ligase_Images")
IMG_SIZE = 800          # px × px high-res output
LINE_WIDTH = 2.0        # thicker bonds

os.makedirs(OUT_DIR, exist_ok=True)

print("🔍 Scanning for SDF files in:", LIGASE_DIR)
sdf_files = [f for f in os.listdir(LIGASE_DIR) if f.endswith(".sdf")]

if not sdf_files:
    print("❌ No SDF files found in Ligases/")
    exit()

print(f"📦 Found {len(sdf_files)} ligase SDFs\n")

generated = 0
failed = 0

for sdf in sdf_files:
    sdf_path = os.path.join(LIGASE_DIR, sdf)
    png_name = sdf.replace(".sdf", ".png")
    png_path = os.path.join(OUT_DIR, png_name)

    try:
        suppl = Chem.SDMolSupplier(sdf_path)
        mol = next(iter(suppl), None)

        if mol is None:
            print(f"⚠️ Could not parse {sdf}, skipping…")
            failed += 1
            continue

        # Generate 2D coordinates if absent
        rdDepictor.Compute2DCoords(mol)

        # High-res RDKit 2D drawer
        drawer = rdMolDraw2D.MolDraw2DCairo(IMG_SIZE, IMG_SIZE)
        drawer.drawOptions().bondLineWidth = LINE_WIDTH
        drawer.drawOptions().addAtomIndices = False
        drawer.drawOptions().addStereoAnnotation = False

        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()

        with open(png_path, "wb") as f:
            f.write(drawer.GetDrawingText())

        print(f"🖼️ Generated: {png_name}")
        generated += 1

    except Exception as e:
        print(f"❌ Error rendering {sdf}: {e}")
        failed += 1

print("\n=============================")
print("🎉 IMAGE GENERATION COMPLETE")
print("=============================")
print(f"🟢 Successfully rendered: {generated}")
print(f"🔴 Failed: {failed}")
print(f"📁 Output folder: {OUT_DIR}")
