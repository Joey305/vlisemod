# DRUGapp.py
from flask import Blueprint, request, jsonify, render_template, current_app
import os
import warnings
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64
import re
import threading
from contextlib import contextmanager

# ---------- Blueprint ----------
# Name it 'dp' so endpoints are 'dp.*' (matches how you register it in app.py)
dp = Blueprint("dp", __name__, template_folder="templates")

# ---------- Env / warnings ----------
os.environ["MKL_SERVICE_FORCE_INTEL"] = "1"
os.environ["MKL_DISABLE_FAST_MM"] = "1"
warnings.filterwarnings("ignore")

# ---------- Concurrency controls ----------
GEN_MAX_CONCURRENT = int(os.environ.get("GEN_MAX_CONCURRENT", 1))     # 1 == serialize on single GPU
GEN_ACQUIRE_TIMEOUT = float(os.environ.get("GEN_ACQUIRE_TIMEOUT", 3)) # seconds
_GEN_SEM = threading.Semaphore(GEN_MAX_CONCURRENT)

@contextmanager
def acquire_gen_slot(timeout=GEN_ACQUIRE_TIMEOUT):
    ok = _GEN_SEM.acquire(timeout=timeout)
    try:
        yield ok
    finally:
        if ok:
            _GEN_SEM.release()

DEFAULT_MODEL_ID = "meta-llama/Llama-3.1-8B-Instruct"


def _env_flag(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def _hf_token():
    return (
        os.environ.get("HF_TOKEN")
        or os.environ.get("HUGGINGFACE_HUB_TOKEN")
        or os.environ.get("HUGGING_FACE_HUB_TOKEN")
    )


def _load_tokenizer(model_id: str):
    kwargs = {}
    token = _hf_token()
    if token:
        kwargs["token"] = token
    return AutoTokenizer.from_pretrained(model_id, **kwargs)


def _load_model(model_id: str):
    kwargs = {"device_map": "auto"}
    token = _hf_token()
    if token:
        kwargs["token"] = token

    try:
        return AutoModelForCausalLM.from_pretrained(
            model_id,
            dtype=torch.float16,
            **kwargs,
        )
    except TypeError:
        return AutoModelForCausalLM.from_pretrained(
            model_id,
            torch_dtype=torch.float16,
            **kwargs,
        )


def configure_local_llm(app):
    model_id = os.environ.get("MODEL_ID") or os.environ.get("LLM_MODEL_ID") or DEFAULT_MODEL_ID
    enabled = _env_flag("ENABLE_LOCAL_LLM", default=False)

    app.config["LOCAL_LLM_ENABLED"] = enabled
    app.config["LOCAL_LLM_MODEL_ID"] = model_id
    app.config["LLAMA_MODEL"] = None
    app.config["LLAMA_TOKENIZER"] = None
    app.config["LOCAL_LLM_ERROR"] = None

    if not enabled:
        app.logger.info("Local LLM disabled. Set ENABLE_LOCAL_LLM=true to enable Drug GPT.")
        return

    app.logger.info("Starting local LLM load for %s", model_id)
    try:
        tokenizer = _load_tokenizer(model_id)
        if tokenizer.pad_token is None and tokenizer.eos_token is not None:
            tokenizer.pad_token = tokenizer.eos_token

        model = _load_model(model_id)
        model.eval()
        torch.set_grad_enabled(False)

        app.config["LLAMA_TOKENIZER"] = tokenizer
        app.config["LLAMA_MODEL"] = model
        app.logger.info("Local LLM load completed for %s", model_id)
    except Exception as exc:
        app.config["LOCAL_LLM_ERROR"] = str(exc)
        app.logger.exception("Local LLM load failed for %s", model_id)

def get_model():
    return current_app.config["LLAMA_MODEL"]

def get_tokenizer():
    return current_app.config["LLAMA_TOKENIZER"]


def local_llm_ready() -> bool:
    return bool(
        current_app.config.get("LOCAL_LLM_ENABLED")
        and current_app.config.get("LLAMA_MODEL") is not None
        and current_app.config.get("LLAMA_TOKENIZER") is not None
    )


def local_llm_error_message() -> str:
    if not current_app.config.get("LOCAL_LLM_ENABLED"):
        return "Local LLM is disabled. Set ENABLE_LOCAL_LLM=true to enable Drug GPT."

    error = current_app.config.get("LOCAL_LLM_ERROR")
    if error:
        return f"Local LLM is unavailable right now: {error}"

    return "Local LLM is not ready yet. Please try again shortly."


def local_llm_unavailable_response(status_code=503):
    return jsonify(
        {
            "error": local_llm_error_message(),
            "local_llm_enabled": bool(current_app.config.get("LOCAL_LLM_ENABLED")),
            "model_id": current_app.config.get("LOCAL_LLM_MODEL_ID"),
        }
    ), status_code

# ---------- Helpers ----------
def extract_assistant_reply(generated_text: str) -> str:
    """Extract only the assistant's final response from a chat-style decode."""
    m = re.search(r'assistant\s+', generated_text)
    return generated_text[m.end():].strip() if m else generated_text.strip()

def format_chat_prompt(user_input: str) -> str:
    tokenizer = get_tokenizer()  # Retrieve tokenizer from app config
    return tokenizer.apply_chat_template(
        [
            {"role": "system", "content": "You are a concise medical assistant. Answer clearly and briefly."},
            {"role": "user", "content": user_input}
        ],
        add_generation_prompt=True,
        tokenize=False
    )

def query_pubchem(compound_name):
    try:
        results = pcp.get_compounds(compound_name, 'name')
        return results[0].to_dict() if results else None
    except Exception as e:
        print(f"PubChem error: {e}")
        return None

def visualize_molecule(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        img = Draw.MolToImage(mol)
        buf = BytesIO()
        img.save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception as e:
        print(f"RDKit error: {e}")
        return None

# ---------- Routes ----------
@dp.route("/")
def home():
    return render_template(
        "DRUGindex.html",
        local_llm_enabled=bool(current_app.config.get("LOCAL_LLM_ENABLED")),
        local_llm_ready=local_llm_ready(),
        local_llm_error=current_app.config.get("LOCAL_LLM_ERROR"),
        model_id=current_app.config.get("LOCAL_LLM_MODEL_ID"),
    )

@dp.route("/query", methods=["POST"])
def query_biogpt():
    """Generate a response from the model (GPU-guarded)."""
    if not local_llm_ready():
        return local_llm_unavailable_response()

    try:
        data = request.get_json(silent=True) or {}
        question = (data.get("question") or "").strip()
        if not question:
            return jsonify({"error": "No question provided"}), 400

        # Acquire a GPU generation slot (non-blocking beyond GEN_ACQUIRE_TIMEOUT)
        with acquire_gen_slot() as ok:
            if not ok:
                return jsonify({"error": "Model is busy. Please try again in a few seconds."}), 429

            # Prepare inputs
            tokenizer = get_tokenizer()  # Retrieve tokenizer from app config
            prompt = format_chat_prompt(question)
            inputs = tokenizer(prompt, return_tensors="pt")
            device = "cuda" if torch.cuda.is_available() else "cpu"
            inputs = {k: v.to(device) for k, v in inputs.items()}

            # Generate under inference mode to save VRAM
            with torch.inference_mode():
                model = get_model()  # Retrieve model from app config
                outputs = model.generate(
                    **inputs,
                    max_new_tokens=200,
                    do_sample=True,
                    temperature=0.5,
                    top_k=40,
                    top_p=0.9,
                    repetition_penalty=1.2,
                    return_dict_in_generate=True,
                    output_scores=False,
                )

        raw = tokenizer.decode(outputs.sequences[0], skip_special_tokens=True)
        response = extract_assistant_reply(raw)
        if "(ABSTRACT TRUNCATED AT 250 WORDS)" in response:
            response = response.replace("(ABSTRACT TRUNCATED AT 250 WORDS)", "").strip()

        # Optional enrichment (CPU-side)
        compound = query_pubchem(question)
        if compound:
            response += (
                f"\n\nCompound Info:"
                f"\nName: {compound.get('IUPACName', 'N/A')}"
                f"\nMolecular Formula: {compound.get('MolecularFormula', 'N/A')}"
                f"\nMolecular Weight: {compound.get('MolecularWeight', 'N/A')}\n"
            )
            smiles = compound.get("CanonicalSMILES", "")
            if smiles:
                img64 = visualize_molecule(smiles)
                if img64:
                    response += (
                        "\nMolecular Structure:\n"
                        f"<img src='data:image/png;base64,{img64}' alt='Molecular Structure' />"
                    )

        return jsonify({"response": response})

    except torch.cuda.OutOfMemoryError:
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        return jsonify({"error": "GPU out of memory during generation. Please retry or shorten the question."}), 503
    except Exception as e:
        # Log and return safe error
        print("Error in /drugapp/query:", repr(e))
        return jsonify({"error": "Unexpected error during generation."}), 500
