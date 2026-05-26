import json
import re
import multiprocessing
from itertools import combinations, islice

def format_comparative_reasoning(better, worse):
    lines = []

    lines.append(f"{better['ligand']} scores higher overall ({better['score']}/4) compared to {worse['ligand']} ({worse['score']}/4).")

    if "linker-friendly functional groups" in better["reasoning"]:
        better_fg = re.search(r"Multiple exposed linker-friendly functional groups: (.+?)\.", better["reasoning"])
        if better_fg:
            lines.append(f"{better['ligand']} has more exposed linker-friendly functional groups ({better_fg.group(1)}).")

    if better["avg_contacts"] > worse["avg_contacts"]:
        lines.append(f"It also shows more average high-quality binding contacts ({better['avg_contacts']:.1f} vs {worse['avg_contacts']:.1f}).")

    if better["interaction_density"] > worse["interaction_density"]:
        lines.append(f"Its interaction density is higher ({better['interaction_density']:.4f} vs {worse['interaction_density']:.4f}), indicating tighter binding for its size.")

    if "molecular_weight" in better and "molecular_weight" in worse:
        if better["molecular_weight"] < worse["molecular_weight"]:
            lines.append(f"{better['ligand']} is also smaller in molecular weight ({better['molecular_weight']:.2f} Da vs {worse['molecular_weight']:.2f} Da).")

    return " ".join(lines)

def generate_pairwise_batch(pairs_batch, successful_ligands, score_threshold=0.5, min_diff_contacts=0.5, min_diff_density=0.001):
    contrastive_data = []

    for i, j in pairs_batch:
        ligand_a = successful_ligands[i]
        ligand_b = successful_ligands[j]

        score_diff = ligand_a["score"] - ligand_b["score"]
        if abs(score_diff) >= score_threshold:
            better, worse = (ligand_a, ligand_b) if score_diff > 0 else (ligand_b, ligand_a)
        else:
            contact_diff = ligand_a.get("avg_contacts", 0) - ligand_b.get("avg_contacts", 0)
            if abs(contact_diff) >= min_diff_contacts:
                better, worse = (ligand_a, ligand_b) if contact_diff > 0 else (ligand_b, ligand_a)
            else:
                density_diff = ligand_a.get("interaction_density", 0) - ligand_b.get("interaction_density", 0)
                if abs(density_diff) >= min_diff_density:
                    better, worse = (ligand_a, ligand_b) if density_diff > 0 else (ligand_b, ligand_a)
                else:
                    continue

        contrastive_data.append({
            "input": f"Which is the better PROTAC warhead candidate: {better['ligand']} or {worse['ligand']}? Why?",
            "output": f"{better['ligand']} is the better candidate.\n\n{format_comparative_reasoning(better, worse)}",
            "metadata": {
                "better": better["ligand"],
                "worse": worse["ligand"],
                "score_diff": round(score_diff, 2),
                "contacts_diff": round(better["avg_contacts"] - worse["avg_contacts"], 2),
                "density_diff": round(better["interaction_density"] - worse["interaction_density"], 4)
            }
        })

    return contrastive_data

def chunked_iterable(iterable, size):
    it = iter(iterable)
    return iter(lambda: list(islice(it, size)), [])

if __name__ == "__main__":
    with open("training_data2.json", "r") as f:
        raw_data = json.load(f)

    successful_ligands = []
    seen = set()

    for entry in raw_data:
        if "Adaptability Score for" in entry["output"]:
            match = re.search(r"Score for (\w+): ([\d\.]+)/4", entry["output"])
            if match:
                ligand = match.group(1)
                if ligand in seen:
                    continue
                seen.add(ligand)

                score = float(match.group(2))
                try:
                    avg = float(re.search(r"([\d\.]+) average high-quality contacts", entry["output"]).group(1))
                    density = float(re.search(r"Interaction Density: ([\d\.]+)", entry["output"]).group(1))
                    mw_match = re.search(r"molecular weight \(([\d\.]+)", entry["output"])
                    mw = float(mw_match.group(1)) if mw_match else 0.0
                except:
                    avg, density, mw = 0.0, 0.0, 0.0

                successful_ligands.append({
                    "ligand": ligand,
                    "score": score,
                    "avg_contacts": avg,
                    "interaction_density": density,
                    "molecular_weight": mw,
                    "reasoning": entry["output"]
                })

    pair_indices = list(combinations(range(len(successful_ligands)), 2))
    num_workers = multiprocessing.cpu_count()
    chunks = list(chunked_iterable(pair_indices, len(pair_indices) // num_workers + 1))

    with multiprocessing.Pool(num_workers) as pool:
        results = pool.starmap(generate_pairwise_batch, [(chunk, successful_ligands) for chunk in chunks])

    contrastive_pairs = [pair for sublist in results for pair in sublist]

    with open("training_data_contrastive.json", "w") as f:
        json.dump(contrastive_pairs, f, indent=4)

    print(f"✅ Saved {len(contrastive_pairs)} contrastive training prompts to training_data_contrastive.json")
