import os
import torch
from datasets import load_dataset, load_from_disk
from transformers import (
    AutoTokenizer,
    AutoModelForCausalLM,
    TrainingArguments,
    Trainer,
    BitsAndBytesConfig
)
from peft import prepare_model_for_kbit_training, LoraConfig, get_peft_model

# === Model ID ===
model_id = "meta-llama/Llama-3.1-8B-Instruct"

# === File Paths ===
input_jsonl = "/mnt/c/Users/joey/Documents/WorkingViralDB/training_data_llama_chat_ready.jsonl"
tokenized_path = "/mnt/d/llama/tokenized_llama_dataset"
output_dir = "/mnt/d/llama/llama3_finetuned"
log_dir = os.path.join(output_dir, "logs")

# === Tokenizer ===
tokenizer = AutoTokenizer.from_pretrained(model_id)
tokenizer.pad_token = tokenizer.eos_token

# === Load or Tokenize Dataset ===
if os.path.exists(tokenized_path):
    print(f"✅ Loading tokenized dataset from {tokenized_path}...")
    tokenized_dataset = load_from_disk(tokenized_path)
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

    tokenized_dataset = dataset.map(preprocess, batched=False)
    os.makedirs(tokenized_path, exist_ok=True)
    tokenized_dataset.save_to_disk(tokenized_path)
    print(f"✅ Tokenization complete and saved to {tokenized_path}.")

# === Quantization Setup ===
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_compute_dtype=torch.float16,
    bnb_4bit_use_double_quant=True,
    bnb_4bit_quant_type="nf4"
)

# === Load Model + Prepare for LoRA ===
model = AutoModelForCausalLM.from_pretrained(
    model_id,
    quantization_config=bnb_config,
    torch_dtype=torch.float16,
    device_map="auto"
)
model = prepare_model_for_kbit_training(model)

# === Apply LoRA Adapter ===
peft_config = LoraConfig(
    r=8,
    lora_alpha=16,
    target_modules=["q_proj", "k_proj", "v_proj", "o_proj"],
    lora_dropout=0.05,
    bias="none",
    task_type="CAUSAL_LM"
)
model = get_peft_model(model, peft_config)

# === Training Args ===
training_args = TrainingArguments(
    output_dir=output_dir,
    per_device_train_batch_size=1,
    gradient_accumulation_steps=4,
    num_train_epochs=3,
    learning_rate=2e-4,
    fp16=True,
    logging_dir=log_dir,
    logging_steps=10,
    save_strategy="epoch",
    save_total_limit=2,
    report_to="none",
    resume_from_checkpoint=True  # <=== Enables automatic resume if interrupted
)

# === Trainer Setup ===
trainer = Trainer(
    model=model,
    train_dataset=tokenized_dataset,
    args=training_args,
    tokenizer=tokenizer
)

# === Start Training ===
trainer.train()
