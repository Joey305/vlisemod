import json
import re
import random

# -------------------------------
# Helper functions for metadata extraction
# -------------------------------
def extract_metadata(entry):
    """Extract ligand, score, avg_contacts, interaction_density, and reasoning from a training_data2.json entry."""
    metadata = {}
    m = re.search(r"Score for (\w+): ([\d\.]+)/4", entry["output"])
    if not m:
        return None
    metadata["ligand"] = m.group(1)
    metadata["score"] = float(m.group(2))
    try:
        avg = float(re.search(r"([\d\.]+) average high-quality contacts", entry["output"]).group(1))
        density = float(re.search(r"Interaction Density: ([\d\.]+)", entry["output"]).group(1))
    except Exception:
        avg, density = 0.0, 0.0
    metadata["avg_contacts"] = avg
    metadata["interaction_density"] = density
    # Optionally extract molecular weight if mentioned
    try:
        mw_match = re.search(r"molecular weight \(([\d\.]+)", entry["output"])
        metadata["molecular_weight"] = float(mw_match.group(1)) if mw_match else None
    except Exception:
        metadata["molecular_weight"] = None
    metadata["reasoning"] = entry["output"]
    return metadata

# -------------------------------
# 1. Instruction-Tuned Single-Turn QA
# -------------------------------
def generate_instruction_qa(metadata):
    ligand = metadata["ligand"]
    # Create a question about linker-friendly groups
    q = f"What linker-friendly functional groups are present in {ligand}?"
    # Try to extract the linker group information from the reasoning text
    fg_match = re.search(r"exposed linker-friendly functional groups: ([^\.]+)\.", metadata["reasoning"], re.IGNORECASE)
    a = f"{ligand} has "
    if fg_match:
        a += fg_match.group(1) + "."
    else:
        a += "no reported linker-friendly functional groups."
    return {"type": "instruction_qa", "input": q, "output": a}

# -------------------------------
# 2. "Explain This Output" Rationalization
# -------------------------------
def generate_explain_output(metadata):
    ligand = metadata["ligand"]
    q = (f"Interpret the following PROTAC metrics for {ligand}: "
         f"Score = {metadata['score']}/4, "
         f"Average Contacts = {metadata['avg_contacts']}, "
         f"Interaction Density = {metadata['interaction_density']}. "
         "What do these numbers suggest about its suitability as a PROTAC warhead?")
    a = (f"For {ligand}, a higher score combined with more average high-quality contacts and a higher interaction density "
         "indicates strong binding efficiency and favorable structural features for PROTAC design. "
         "This suggests that it likely has multiple exposed functional groups and robust interaction patterns, making it a promising candidate.")
    return {"type": "explain_output", "input": q, "output": a}

# -------------------------------
# 4. Ranking Task / Preference Modeling
# -------------------------------
def generate_ranking_example(metadata_list):
    # Select three distinct ligands at random
    if len(metadata_list) < 3:
        return None
    candidates = random.sample(metadata_list, 3)
    ligand_names = [m["ligand"] for m in candidates]
    q = f"Rank the following PROTAC warhead candidates in order of overall suitability and explain your reasoning: {', '.join(ligand_names)}."
    
    # Sort by score then by avg_contacts (as a secondary criterion)
    ranked = sorted(candidates, key=lambda m: (m["score"], m["avg_contacts"]), reverse=True)
    ranking_str = ", ".join([f"{m['ligand']} (Score: {m['score']}/4)" for m in ranked])
    a = f"Based on the provided metrics, the ranking is: {ranking_str}. " \
        "This ordering is determined by the overall score, which reflects the combined impact of binding interactions and structural features."
    return {"type": "ranking", "input": q, "output": a}

# -------------------------------
# 5. Chain-of-Thought (CoT) Expansion
# -------------------------------
def generate_chain_of_thought(metadata):
    ligand = metadata["ligand"]
    q = f"Explain step-by-step why {ligand} is a good PROTAC candidate."
    a = (
        f"Step 1: {ligand} has a high overall score ({metadata['score']}/4), indicating strong binding interactions.\n"
        f"Step 2: It forms an average of {metadata['avg_contacts']} high-quality contacts, which shows robust molecular interactions.\n"
        f"Step 3: Its interaction density of {metadata['interaction_density']:.4f} suggests that it binds efficiently relative to its size.\n"
        "Therefore, these factors combined indicate that it is a promising candidate for PROTAC design."
    )
    return {"type": "chain_of_thought", "input": q, "output": a}

# -------------------------------
# 8. Synthetic Negative Reasoning
# -------------------------------
def generate_negative_reasoning(metadata):
    ligand = metadata["ligand"]
    # Use a threshold; if score is low (< 2), consider it a poor candidate
    if metadata["score"] >= 2:
        return None
    q = f"Is {ligand} a suitable candidate for PROTAC development? Explain why or why not."
    a = (f"No, {ligand} is not a suitable candidate for PROTAC development. "
         f"It only achieves a score of {metadata['score']}/4, which suggests it has insufficient high-quality binding interactions and "
         "lacks adequate exposure of linker-friendly functional groups. These shortcomings mean it may not be able to form a stable ternary complex.")
    return {"type": "negative_reasoning", "input": q, "output": a}

# -------------------------------
# 9. Glossary / Concept Definition
# -------------------------------
def generate_glossary_examples():
    examples = []
    glossary = {
        "Linker-friendly functional group": (
            "A chemical group (such as hydroxyl, amine, carboxylic acid, or thiol) that is accessible on the surface of a ligand, "
            "making it a suitable attachment point for a linker in PROTAC design."
        ),
        "Interaction Density": (
            "A measure of how many high-quality binding contacts a ligand makes per Dalton of molecular weight, "
            "indicating how efficiently it binds relative to its size."
        ),
        "Molecular Weight Threshold": (
            "An optimal range for PROTAC ligands, typically lower than 500 Da, ensures that the ligand is sufficiently small "
            "to allow for effective linker and binder design."
        )
    }
    for term, definition in glossary.items():
        q = f"What is meant by '{term}'?"
        a = definition
        examples.append({"type": "glossary", "input": q, "output": a})
    return examples

# -------------------------------
# Main script
# -------------------------------
if __name__ == "__main__":
    with open("training_data2.json", "r") as f:
        raw_data = json.load(f)

    metadata_list = []
    for entry in raw_data:
        meta = extract_metadata(entry)
        if meta:
            metadata_list.append(meta)

    extended_examples = []

    # 1. Instruction-Tuned Single-Turn QA (for each ligand)
    for meta in metadata_list:
        extended_examples.append(generate_instruction_qa(meta))

    # 2. "Explain This Output" Rationalization (for each ligand)
    for meta in metadata_list:
        extended_examples.append(generate_explain_output(meta))

    # 4. Ranking Task / Preference Modeling (sample a few groups)
    num_ranking = min(5, len(metadata_list) // 3)
    for _ in range(num_ranking):
        ranking_example = generate_ranking_example(metadata_list)
        if ranking_example:
            extended_examples.append(ranking_example)

    # 5. Chain-of-Thought Expansion (for each ligand)
    for meta in metadata_list:
        extended_examples.append(generate_chain_of_thought(meta))

    # 8. Synthetic Negative Reasoning (for ligands with low score)
    for meta in metadata_list:
        neg_example = generate_negative_reasoning(meta)
        if neg_example:
            extended_examples.append(neg_example)

    # 9. Glossary / Concept Definitions
    extended_examples.extend(generate_glossary_examples())

    # Shuffle the extended examples for variety
    random.shuffle(extended_examples)

    with open("training_data_extended.json", "w") as f:
        json.dump(extended_examples, f, indent=4)

    print(f"✅ Saved {len(extended_examples)} extended training examples to extended_training_data.json")
