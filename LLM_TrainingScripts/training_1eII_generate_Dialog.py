import json
import re
import multiprocessing
from itertools import islice

# Challenge questions to simulate skepticism or tradeoffs
CHALLENGE_TEMPLATES = [
    "But what if the other ligand has better pharmacokinetics?",
    "What about synthetic accessibility for the other molecule?",
    "Could the larger molecule offer better selectivity in vivo?",
    "Doesn't a lower molecular weight usually help with oral bioavailability?",
    "Could the other ligand be easier to modify for linker attachment?",
    "What if the weaker one is more stable in plasma?",
    "Could the size of the better candidate be a disadvantage?",
    "But what if the other ligand is more specific to the target?",
]

def format_followup(better, worse, challenge):
    return {
        "role": "assistant",
        "content": (
            f"That's a good point. While {worse} may offer that advantage, "
            f"{better} still remains the better candidate due to its stronger interaction profile, "
            f"exposed linker groups, and tighter binding footprint."
        )
    }

def generate_dialog_from_pair(pair, challenge_templates):
    dialog = []

    # Extract ligand names and winner
    q_match = re.search(r"Which is the better PROTAC warhead candidate: (\w+) or (\w+)", pair["input"])
    a_match = re.search(r"^(\w+) is the better candidate", pair["output"])
    if not q_match or not a_match:
        return None

    ligand_a, ligand_b = q_match.groups()
    better = a_match.group(1)
    worse = ligand_b if better == ligand_a else ligand_a

    # Base Q&A
    dialog.append({"role": "user", "content": pair["input"]})
    dialog.append({"role": "assistant", "content": pair["output"]})

    # Add one challenge Q+A
    challenge = challenge_templates[hash(better + worse) % len(challenge_templates)]
    dialog.append({"role": "user", "content": challenge})
    dialog.append(format_followup(better, worse, challenge))

    return {"dialog": dialog}

def chunked_iterable(iterable, size):
    it = iter(iterable)
    return iter(lambda: list(islice(it, size)), [])

def process_chunk(chunk):
    return [dialog for dialog in (generate_dialog_from_pair(pair, CHALLENGE_TEMPLATES) for pair in chunk) if dialog]

if __name__ == "__main__":
    with open("training_data_contrastive.json", "r") as f:
        contrastive_pairs = json.load(f)

    num_workers = multiprocessing.cpu_count()
    chunks = list(chunked_iterable(contrastive_pairs, len(contrastive_pairs) // num_workers + 1))

    with multiprocessing.Pool(num_workers) as pool:
        results = pool.map(process_chunk, chunks)

    dialogs = [dialog for batch in results for dialog in batch]

    with open("training_data_contrastive_dialog.json", "w") as f:
        json.dump(dialogs, f, indent=4)

    print(f"💬 Saved {len(dialogs)} multi-turn dialog examples to training_data_contrastive_dialog.json")
