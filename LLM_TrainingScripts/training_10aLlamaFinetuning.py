import os
import time
from datasets import load_dataset
from transformers import AutoTokenizer

# === Timer Start ===
start_time = time.time()

# === Model ID ===
model_id = "meta-llama/Llama-3.1-8B-Instruct"

# === File Paths ===
input_jsonl = "/mnt/d/llama/training_data_llama_chat_ready.jsonl"
tokenized_path = "/mnt/d/llama/tokenized_llama_dataset"

# === Tokenizer ===
print("🔑 Loading tokenizer...")
tokenizer = AutoTokenizer.from_pretrained(model_id)
tokenizer.pad_token = tokenizer.eos_token

# === Load and Tokenize ===
if os.path.exists(tokenized_path) and os.listdir(tokenized_path):
    print(f"✅ Tokenized dataset already exists at: {tokenized_path}")
else:
    print("🌀 Tokenizing raw dataset...")
    dataset = load_dataset("json", data_files=input_jsonl, split="train")

    def preprocess(example):
        prompt = tokenizer.apply_chat_template(example["messages"], tokenize=False)
        encoded = tokenizer(prompt, padding="max_length", truncation=True, max_length=2048)
        return {
            "input_ids": encoded["input_ids"],
            "attention_mask": encoded["attention_mask"],
            "labels": encoded["input_ids"]
        }

    tokenized_dataset = dataset.map(
        preprocess,
        batched=True,
        batch_size=512,
        num_proc=os.cpu_count()
    )

    os.makedirs(tokenized_path, exist_ok=True)
    tokenized_dataset.save_to_disk(tokenized_path)
    print(f"✅ Tokenization complete and saved to: {tokenized_path}")

# === Timer End ===
elapsed = time.time() - start_time
print(f"⏱️ Total time: {elapsed / 60:.2f} minutes")
