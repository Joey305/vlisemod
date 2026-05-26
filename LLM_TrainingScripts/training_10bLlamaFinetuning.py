import os
import torch
import numpy as np
from datasets import load_from_disk
from transformers import (
    AutoTokenizer,
    AutoModelForCausalLM,
    TrainingArguments,
    Trainer,
    BitsAndBytesConfig
)
from peft import prepare_model_for_kbit_training, LoraConfig, get_peft_model

# === Paths & Model ===
model_id = "meta-llama/Llama-3.1-8B-Instruct"
tokenized_path = "/mnt/d/llama/tokenized_llama_dataset"
output_dir = "/mnt/d/llama/llama3_finetuned"
log_dir = os.path.join(output_dir, "logs")

# === Load Dataset ===
print(f"✅ Loading tokenized dataset from {tokenized_path}...")
tokenized_dataset = load_from_disk(tokenized_path)

# === Tokenizer ===
tokenizer = AutoTokenizer.from_pretrained(model_id)
tokenizer.pad_token = tokenizer.eos_token

# === Quantization Setup ===
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_compute_dtype=torch.float16,
    bnb_4bit_use_double_quant=True,
    bnb_4bit_quant_type="nf4"
)

# === Load Base Model & Apply LoRA ===
model = AutoModelForCausalLM.from_pretrained(
    model_id,
    quantization_config=bnb_config,
    torch_dtype=torch.float16,
    device_map="auto"
)
model.config.use_cache = False  # Required for gradient checkpointing
model = prepare_model_for_kbit_training(model)

peft_config = LoraConfig(
    r=8,
    lora_alpha=16,
    target_modules=["q_proj", "k_proj", "v_proj", "o_proj"],
    lora_dropout=0.05,
    bias="none",
    task_type="CAUSAL_LM"
)
model = get_peft_model(model, peft_config)

# === Fix UnpicklingError from RNG (torch 2.2+) ===
torch.serialization.add_safe_globals({
    np.dtype: "numpy.dtype",
    np._core.multiarray._reconstruct: "numpy._core.multiarray._reconstruct",
    np.ndarray: "numpy.ndarray",
    np.dtypes.UInt32DType: "numpy.dtypes.UInt32DType"
})

# === Training Arguments ===
training_args = TrainingArguments(
    output_dir=output_dir,
    per_device_train_batch_size=1,
    gradient_accumulation_steps=4,
    num_train_epochs=3,
    learning_rate=2e-4,
    fp16=True,
    gradient_checkpointing=True,
    logging_dir=log_dir,
    logging_steps=10,
    save_strategy="steps",
    save_steps=100,
    save_total_limit=2,
    report_to="none"
)

# === Setup Trainer ===
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_dataset,
    tokenizer=tokenizer
)

# === Auto Resume Logic ===
checkpoints = [
    os.path.join(output_dir, d)
    for d in os.listdir(output_dir)
    if d.startswith("checkpoint-") and os.path.isdir(os.path.join(output_dir, d))
]

if checkpoints:
    latest_checkpoint = max(checkpoints, key=lambda x: int(x.split("-")[-1]))
    print(f"🔁 Resuming from latest checkpoint: {latest_checkpoint}")
    try:
        trainer.train(resume_from_checkpoint=latest_checkpoint)
    except Exception as e:
        print(f"⚠️ Failed to load RNG state from checkpoint: {e}")
        print("🧯 Retrying without loading RNG...")
        trainer.train(resume_from_checkpoint=False)
else:
    print("🆕 Starting training from scratch.")
    trainer.train()
