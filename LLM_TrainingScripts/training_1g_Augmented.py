import json
import re
from tqdm import tqdm

# Load the contrastive dataset
with open("training_data_contrastive.json", "r") as f:
    contrastive_data = json.load(f)

# Generate new reasoning examples from top 4 reasoning types

def generate_reasoning_prompts(pair):
    prompts = []

    q_match = re.search(r"Which is the better PROTAC warhead candidate: (\w+) or (\w+)", pair["input"])
    a_match = re.search(r"^(\w+) is the better candidate", pair["output"])
    
    if not q_match or not a_match:
        return []

    ligand_a, ligand_b = q_match.groups()
    better = a_match.group(1)
    worse = ligand_b if better == ligand_a else ligand_a

    reasoning = pair["output"].split("\n\n")[-1]

    # 1. Reverse justification
    prompts.append({
        "input": f"Why not use {worse} instead?",
        "output": f"While {worse} may have some advantages, {better} is superior overall. {reasoning}"
    })

    # 2. Design prompt
    prompts.append({
        "input": f"What would improve {better} as a PROTAC warhead?",
        "output": (
            f"To improve {better}, one could increase the number of solvent-exposed linker-friendly groups or optimize its binding affinity. "
            f"Currently, its strength lies in: {reasoning}"
        )
    })

    # 3. Explain like I'm 5
    prompts.append({
        "input": f"Why is {better} good? (Explain like I'm 5)",
        "output": (
            f"{better} sticks really well to the target and has parts we can grab to make special medicines. {worse} doesn't do this as nicely."
        )
    })

    # 4. Critique this reasoning
    prompts.append({
        "input": f"Is this a good explanation? '{pair['output']}'",
        "output": (
            f"It's a good start and explains why {better} is better, but could be improved by comparing molecular weight or interaction density more clearly."
        )
    })

    # 5. Role-Reversal: Improve the weaker one
    prompts.append({
        "input": f"What would make {worse} a better PROTAC warhead than {better}?",
        "output": (
            f"If {worse} had more solvent-exposed functional groups or improved binding affinity, "
            f"it might outperform {better}. Currently, it lacks in comparison."
        )
    })

    # 6. One-sentence summary
    prompts.append({
        "input": f"Summarize why {better} is better than {worse} in one sentence.",
        "output": f"{better} has stronger binding, more exposed linker groups, and better interaction efficiency than {worse}."
    })

    # 7. True or False quiz
    prompts.append({
        "input": f"True or False: {worse} is the better PROTAC candidate compared to {better}.",
        "output": f"False. {better} is superior because {reasoning}"
    })

    return prompts


# Process all pairs
all_reasoning = []
for entry in tqdm(contrastive_data):
    all_reasoning.extend(generate_reasoning_prompts(entry))

# Save
with open("training_data_reasoning_prompts.json", "w") as f:
    json.dump(all_reasoning, f, indent=4)

print(f"✅ Saved {len(all_reasoning)} reasoning-style prompts to training_data_reasoning_prompts.json")
