"""
direct_prompt_baseline.py
=========================
Direct prompting baseline for wheat gene set reconstruction.

Asks the LLM directly: "Given this GO term, which Triticum aestivum genes
belong in this gene set?" — no retrieval, no grading, no documents.

Serves as an evaluation baseline to measure how much the RAG pipeline
adds over zero-shot LLM knowledge.

Outputs (mirroring the maker pipeline structure):
  out/geneset_generations_direct/{llm_name}/{gene_set_name}.json
  out/geneset_generations_direct/{llm_name}/{gene_set_name}_raw.txt
  out/processed_genesets_direct_{llm_name}.txt

Usage:
    # Ollama (default)
    python direct_prompt_baseline.py --llm deepseek-r1:8b

    # Groq
    python direct_prompt_baseline.py --llm llama-3.1-8b-instant --backend groq
"""

import os
import json
import argparse
import re
import time

from langchain_core.messages import HumanMessage, SystemMessage

from utils import phenotype_json_reader

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

OUTPUT_BASE    = "out/geneset_generations_direct"
PROCESSED_BASE = "out/processed_genesets_direct"
DEFAULT_INPUT  = "out/go_terms.json"
RETRY_DELAY    = 2.0
MAX_RETRIES    = 3

SYSTEM_PROMPT = """You are a wheat genomics expert with deep knowledge of Triticum aestivum gene functions.

You will be given a GO (Gene Ontology) term name and its definition. Your task is to list all
Triticum aestivum genes you know that are associated with this GO term.

Return ONLY a JSON array of gene objects. Each object must have exactly these fields:
  "Gene"             — the TraesCS stable ID (e.g. TraesCS4B03G0844400) if known.
                       If you do not know the TraesCS ID, use the functional name from literature (e.g. TaHSP70).
                       Always prefer the TraesCS ID over the functional name.
  "Source Reference" — a brief note from your knowledge supporting the association
  "PMID"             — leave as empty string ""
  "Journal"          — leave as empty string ""

If a gene has both a functional name and a TraesCS ID, use the TraesCS ID in the Gene field
and put the functional name in the Source Reference field.

If you do not know of any wheat genes for this GO term, return an empty array: []
Return ONLY the JSON array. No explanation, no markdown fences.
"""


# ---------------------------------------------------------------------------
# LLM factory
# ---------------------------------------------------------------------------

def get_llm(llm_name: str, backend: str):
    """
    Return a LangChain chat model for the given backend.
    Ollama: uses ChatOllama
    Groq:   uses ChatGroq (requires GROQ_API_KEY env var)
    """
    if backend == "groq":
        from langchain_groq import ChatGroq
        return ChatGroq(
            model=llm_name,
            api_key=os.environ.get("GROQ_API_KEY"),
            temperature=0,
            max_tokens=2048,
        )
    else:
        from langchain_ollama import ChatOllama
        return ChatOllama(
            model=llm_name,
            temperature=0,
        )


# ---------------------------------------------------------------------------
# JSON helpers (mirrors rag_pipeline_gene_set_maker_llm.py)
# ---------------------------------------------------------------------------

def repair_json(raw: str) -> list:
    raw = re.sub(r"<think>.*?</think>", "", raw, flags=re.DOTALL).strip()
    match = re.search(r"```(?:json)?\s*([\s\S]*?)```", raw)
    if match:
        raw = match.group(1).strip()
    match = re.search(r"(\[[\s\S]*\])", raw)
    if match:
        raw = match.group(1).strip()
    if raw.count("[") > raw.count("]"):
        raw += "]"
    if raw.count("{") > raw.count("}"):
        raw += "}"
    return json.loads(raw)


def safe_parse(raw: str) -> list:
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        pass
    try:
        return repair_json(raw)
    except Exception:
        return []


def append_safe_save(new_genes: list, path: str) -> int:
    existing = []
    if os.path.exists(path):
        try:
            with open(path) as f:
                existing = json.load(f)
            if not isinstance(existing, list):
                existing = []
        except Exception:
            existing = []

    seen  = {(str(g.get("Gene", "")).strip().lower(), str(g.get("PMID", "")).strip()) for g in existing}
    added = 0
    for g in new_genes:
        key = (str(g.get("Gene", "")).strip().lower(), str(g.get("PMID", "")).strip())
        if key not in seen:
            existing.append(g)
            seen.add(key)
            added += 1

    with open(path, "w") as f:
        json.dump(existing, f, indent=2)
    return added


def load_processed(path: str) -> set:
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def mark_processed(name: str, path: str):
    with open(path, "a") as f:
        f.write(f"{name}\n")


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------

def generate_direct(geneset: dict, llm, backend: str) -> list:
    """
    Ask the LLM directly for wheat genes in a GO term.
    Retries on rate limit / error with backoff.
    """
    safe_name  = geneset["raw_go_name"].replace("/", " or ")
    definition = geneset.get("go_term_definition", "")

    messages = [
        SystemMessage(content=SYSTEM_PROMPT),
        HumanMessage(content=(
            f"GO term: {safe_name}\n"
            f"Definition: {definition}\n\n"
            f"List all Triticum aestivum genes associated with this GO term. "
            f"Use TraesCS stable IDs (e.g. TraesCS4B03G0844400) wherever possible. "
            f"Only use functional names (e.g. TaBX1) if you do not know the TraesCS ID."
        )),
    ]

    for attempt in range(MAX_RETRIES):
        try:
            result = llm.invoke(messages)
            genes  = safe_parse(result.content)
            print(f"  [{backend}] {len(genes)} gene(s) returned")
            return genes

        except Exception as e:
            err_str = str(e)
            # Handle Groq rate limit — parse wait time from error message
            wait_match = re.search(r"try again in ([0-9.]+)s", err_str)
            if wait_match:
                wait = float(wait_match.group(1)) + 1.0
                print(f"  [rate limit] waiting {wait:.1f}s (attempt {attempt+1}/{MAX_RETRIES})")
                time.sleep(wait)
            else:
                print(f"  [error] attempt {attempt+1}/{MAX_RETRIES}: {e}")
                time.sleep(RETRY_DELAY)

    return []


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Direct prompting baseline: ask LLM for wheat genes per GO term"
    )
    parser.add_argument("--llm",        default="deepseek-r1:8b",
                        help="LLM model name (default: deepseek-r1:8b)")
    parser.add_argument("--backend",    choices=["ollama", "groq"], default="ollama",
                        help="LLM backend to use (default: ollama)")
    parser.add_argument("--input_file", default=DEFAULT_INPUT,
                        help=f"GO terms JSON file (default: {DEFAULT_INPUT})")
    args = parser.parse_args()

    llm_name = args.llm
    backend  = args.backend

    out_dir        = os.path.join(OUTPUT_BASE, llm_name)
    processed_file = f"{PROCESSED_BASE}_{llm_name}.txt"
    os.makedirs(out_dir, exist_ok=True)

    llm       = get_llm(llm_name, backend)
    genesets  = phenotype_json_reader(args.input_file)
    processed = load_processed(processed_file)

    to_process = [g for g in genesets if g["raw_go_name"].replace("/", " or ") not in processed]

    print(f"Total gene sets : {len(genesets)}")
    print(f"Already done    : {len(processed)}")
    print(f"To process      : {len(to_process)}")
    print(f"LLM             : {llm_name} via {backend}\n")

    for geneset in to_process:
        safe_name    = geneset["raw_go_name"].replace("/", " or ")
        json_outfile = os.path.join(out_dir, f"{safe_name}.json")
        raw_outfile  = os.path.join(out_dir, f"{safe_name}_raw.txt")

        print(f"\nProcessing: {safe_name}")

        if not os.path.exists(json_outfile):
            with open(json_outfile, "w") as f:
                json.dump([], f)

        try:
            genes = generate_direct(geneset, llm, backend)

            with open(raw_outfile, "w") as f:
                f.write(json.dumps(genes, indent=2))

            added = append_safe_save(genes, json_outfile)
            print(f"  {len(genes)} gene(s) found, {added} new → {json_outfile}")
            mark_processed(safe_name, processed_file)

        except Exception as e:
            print(f"  [error] {safe_name}: {e}")
            continue

    print(f"\n=== Done ===")
    print(f"  Outputs  → {out_dir}")
    print(f"  Progress → {processed_file}")


if __name__ == "__main__":
    main()