import json
import multiprocessing
import re
from itertools import islice


def expand_to_multiturn(entry):
    turns = []

    try:
        output_section = entry["output"].split("\n\n")[-1]
        lines = output_section.split(". ")
        input_q = entry["input"]

        turns.append({
            "input": input_q,
            "output": entry["output"]
        })

        for line in lines:
            line = line.strip().rstrip(".")
            if "scores higher overall" in line:
                turns.append({
                    "input": f"Why does {line.split()[0]} score higher overall?",
                    "output": line + "."
                })
            elif "linker-friendly functional groups" in line:
                turns.append({
                    "input": "What does it mean to have more linker-friendly functional groups?",
                    "output": line + "."
                })
            elif "average high-quality binding contacts" in line:
                turns.append({
                    "input": "Why do more high-quality contacts matter for PROTACs?",
                    "output": line + "."
                })
            elif "interaction density" in line:
                turns.append({
                    "input": "What is interaction density and why is it important?",
                    "output": line + "."
                })

    except Exception as e:
        print(f"[!] Skipping malformed entry: {e}")
        return []

    return turns


def chunked_iterable(iterable, size):
    it = iter(iterable)
    return iter(lambda: list(islice(it, size)), [])


def process_batch(batch):
    result = []
    for entry in batch:
        result.extend(expand_to_multiturn(entry))
    return result


if __name__ == "__main__":
    with open("training_data_contrastive.json", "r") as f:
        contrastive_data = json.load(f)

    num_workers = multiprocessing.cpu_count()
    chunks = list(chunked_iterable(contrastive_data, len(contrastive_data) // num_workers + 1))

    with multiprocessing.Pool(num_workers) as pool:
        results = pool.map(process_batch, chunks)

    multi_turn_data = [entry for sublist in results for entry in sublist]

    with open("training_data_contrastive_multiturn.json", "w") as f:
        json.dump(multi_turn_data, f, indent=4)

    print(f"🧠 Saved {len(multi_turn_data)} multi-turn Q&A examples to training_data_contrastive_multiturn.json")
