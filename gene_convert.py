"""
convert_genes.py
----------------
Pre-BLAST gene conversion pipeline.

For each gene extracted by the maker LLM:
  1. Try BioMart chip table lookup (instant)
  2. If not found: tag as 'unconverted_gene' (queued for BLAST)
  3. Write conversion_results.json + lightweight HTML viewer

All TraesCS IDs are upgraded to v2.1 (03G) using v11_to_v21_mapping.json.

# No LLM (same as before)
python gene_convert.py --llm deepseek-r1:8b

# LLM via Ollama (default model: deepseek-r1:8b)
python gene_convert.py --llm deepseek-r1:8b --use-llm --llm-backend ollama

# LLM via Ollama with a specific model
python gene_convert.py --llm deepseek-r1:8b --use-llm --llm-backend ollama --llm-model llama3.2:3b

# LLM via Groq (needs GROQ_API_KEY env var, default model: llama3-8b-8192)
python gene_convert.py --llm deepseek-r1:8b --use-llm --llm-backend groq

# Groq with specific model
python gene_convert.py --llm deepseek-r1:8b --use-llm --llm-backend groq --llm-model mixtral-8x7b-32768

Usage:
    python gene_convert.py --llm deepseek-r1:8b
    python gene_convert.py --llm deepseek-r1:8b --biomart out/merged_biomart_with_v11.tsv
"""

import os
import re
import json
import csv
import argparse
import time
from pathlib import Path
from collections import defaultdict

# ── config ────────────────────────────────────────────────────────────────────
DEFAULT_BIOMART        = "out/merged_biomart_with_v11.tsv"
DEFAULT_GEN_BASE       = "out/merged/geneset_generations"
OUTPUT_DIR             = "out/conversion"
RESULTS_JSON           = "out/conversion/conversion_results.json"
VIEWER_HTML            = "out/conversion/viewer.html"
MAKER_FILE             = "out/merged/generated_gene_sets.txt"
V11_TO_V21_MAPPING     = "v11_to_v21_mapping.json"
FUNPLANTGENES_MAPPING  = "out/funplantgenes_mapping.json"

DEFAULT_LLM_MODEL_GROQ   = "llama3-8b-8192"
DEFAULT_LLM_MODEL_OLLAMA = "deepseek-r1:8b"
LLM_RETRY_DELAY          = 1.0   # seconds between retries on API error

# Tags
TAG_CONVERTED        = "converted"
TAG_UNCONVERTED_GENE = "unconverted_gene"
TAG_NOT_GENE         = "not_gene"
# ─────────────────────────────────────────────────────────────────────────────


# ── v1.1 → v2.1 mapping ──────────────────────────────────────────────────────

def load_v21_mapping(path: str = V11_TO_V21_MAPPING) -> dict:
    """Load v1.1 → v2.1 TraesCS ID mapping."""
    if not Path(path).exists():
        print(f"WARNING: v11→v21 mapping not found: {path} — IDs will not be upgraded to v2.1")
        return {}
    with open(path) as f:
        mapping = json.load(f)
    print(f"v11→v21 mapping loaded: {len(mapping):,} entries")
    return mapping


def load_funplantgenes_mapping(path: str = FUNPLANTGENES_MAPPING) -> dict:
    """
    Load funplantgenes Ta-name → TraesCS v1.1 ID mapping.
    Keys are functional gene names (e.g. 'ALI-1', 'FT-A1').
    Values are v1.1 TraesCS IDs — will be upgraded to v2.1 via to_v21().
    """
    if not Path(path).exists():
        print(f"WARNING: funplantgenes mapping not found: {path} — skipping")
        return {}
    with open(path) as f:
        mapping = json.load(f)
    print(f"funplantgenes mapping loaded: {len(mapping):,} entries")
    return mapping



def to_v21(traescs_id: str, v21_mapping: dict) -> str:
    """
    Upgrade a TraesCS ID to v2.1 if a mapping exists.
    Strips transcript suffix and LC suffix before lookup.
    Returns the v2.1 ID if found, else the cleaned input.
    """
    base = traescs_id.split(".")[0]                         # strip .1 .2 etc
    base = re.sub(r"LC$", "", base, flags=re.I)             # strip LC suffix
    return v21_mapping.get(base, base)


# ── BioMart lookup table ──────────────────────────────────────────────────────

def build_lookup_tables(biomart_file: str) -> dict:
    """
    Build in-memory lookup dicts from the merged BioMart TSV.
    Every identifier type maps to Gene stable ID (TraesCS...).
    """
    if not Path(biomart_file).exists():
        print(f"WARNING: BioMart file not found: {biomart_file}")
        return {}

    lookups = defaultdict(dict)
    col_map = {
        "transcript_id":  "Transcript stable ID",
        "gene_name":      "Gene name",
        "gene_synonym":   "Gene Synonym",
        "ncbi_id":        "NCBI gene (formerly Entrezgene) ID",
        "ncbi_accession": "NCBI gene (formerly Entrezgene) accession",
        "refseq_mrna":    "RefSeq mRNA ID",
        "uniprot_swiss":  "UniProtKB/Swiss-Prot ID",
        "uniprot_trembl": "UniProtKB/TrEMBL ID",
        "uniprot_symbol": "UniProtKB Gene Name symbol",
        "wikigene":       "WikiGene name",
        "interpro_desc":  "InterPro Short Description",
        "protein_domain": "Protein domain description",
        "ncbi_gene_desc": "NCBI gene (formerly Entrezgene) description",
    }

    with open(biomart_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row.get("Gene stable ID", "").strip()
            if not gene_id or not gene_id.startswith("TraesCS"):
                continue
            for key, col in col_map.items():
                val = row.get(col, "").strip()
                if val and val not in ("", "nan"):
                    lookups[key][val.lower()] = gene_id
                    lookups[key][val]          = gene_id

    total = sum(len(v) for v in lookups.values())
    print(f"Lookup tables built: {len(lookups)} types, {total:,} total entries")
    return dict(lookups)


def chip_lookup(name: str, lookups: dict, v21_mapping: dict, funplantgenes: dict) -> tuple[str | None, str]:
    """
    Try all lookup tables. Returns (TraesCS_v2.1_ID, source) or (None, 'not_found').
    All returned IDs are upgraded to v2.1 via v21_mapping.

    Order:
      1. Direct TraesCS ID — upgrade to v2.1
      2. funplantgenes Ta-name lookup (v1.1 → v2.1)
      3. BioMart tables
    """
    name_lower = name.strip().lower()
    name_exact = name.strip()

    # Already in TraesCS format — upgrade to v2.1 and return
    if re.match(r"^TraesCS\w+", name_exact):
        return to_v21(name_exact, v21_mapping), "direct_traescs"

    # funplantgenes Ta-name lookup (case-insensitive)
    fpg_hit = funplantgenes.get(name_exact) or funplantgenes.get(name_lower)
    if fpg_hit:
        return to_v21(fpg_hit, v21_mapping), "funplantgenes"

    # old Traes_ MIPS format and all others — fall through to BioMart tables

    order = [
        ("gene_name",      "biomart_gene_name"),
        ("gene_synonym",   "biomart_synonym"),
        ("transcript_id",  "biomart_transcript"),
        ("refseq_mrna",    "biomart_refseq"),
        ("uniprot_swiss",  "biomart_uniprot_swiss"),
        ("uniprot_symbol", "biomart_uniprot_symbol"),
        ("wikigene",       "biomart_wikigene"),
        ("ncbi_id",        "biomart_ncbi"),
        ("ncbi_accession", "biomart_ncbi_accession"),
        ("uniprot_trembl", "biomart_trembl"),
        ("interpro_desc",  "biomart_interpro"),
        ("protein_domain", "biomart_protein_domain"),
    ]

    for table_key, source_label in order:
        table  = lookups.get(table_key, {})
        result = table.get(name_exact) or table.get(name_lower)
        if result:
            return to_v21(result, v21_mapping), source_label

    return None, "not_found"


# ── LLM normalisation (optional) ─────────────────────────────────────────────

# ── LLM normalisation (optional) ─────────────────────────────────────────────

LLM_CACHE_FILE  = "out/conversion/llm_normalise_cache.json"
LLM_BATCH_SIZE  = 40

# Patterns that are definitely not individual genes — skip LLM entirely
NON_GENE_PATTERNS = [
    r"^\d+(\.\d+)?$",                          # pure numbers
    r"^(wheat|triticum|rice|arabidopsis|maize|barley|sorghum)$",  # organisms
    r"\b(pathway|process|complex|family|domain|motif|class|group|type|clade)\b",
    r"^(high|low|wild|mutant|control|treatment|sample|line|cultivar|variety)$",
    r"^\w{1,2}$",                               # single/double char tokens
]
NON_GENE_RE = re.compile("|".join(NON_GENE_PATTERNS), re.IGNORECASE)


def is_obvious_non_gene(name: str) -> bool:
    return bool(NON_GENE_RE.search(name))


def load_llm_cache() -> dict:
    p = Path(LLM_CACHE_FILE)
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {}


def save_llm_cache(cache: dict):
    Path(LLM_CACHE_FILE).parent.mkdir(parents=True, exist_ok=True)
    with open(LLM_CACHE_FILE, "w") as f:
        json.dump(cache, f, indent=2)


LLM_BATCH_SYSTEM_PROMPT = """You are a wheat genomics expert helping normalise gene names for BioMart lookup.

You will receive a JSON array of gene name strings. For each one:
1. Decide if it is actually a gene (not a pathway, complex, organism, method, or non-gene term).
2. If it is a gene, return the most BioMart-searchable form of the name.

Respond ONLY with a JSON array of objects in the same order as the input, like:
[
  {"is_gene": true, "normalised": "HSP70"},
  {"is_gene": false, "normalised": null},
  ...
]

Rules:
- Strip species prefixes only if unambiguous (TaHSP70 -> HSP70, TaWRKY1 -> WRKY1)
- Keep subgenome suffixes (-A, -B, -D) when part of the name
- If already a TraesCS ID, return is_gene=true and normalised=the input unchanged
- If unsure, return is_gene=false
- Output ONLY the JSON array, no explanation
"""


def _call_groq_batch(names: list[str], model: str) -> list[dict] | None:
    try:
        import groq
        client   = groq.Groq(api_key=os.environ.get("GROQ_API_KEY"))
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": LLM_BATCH_SYSTEM_PROMPT},
                {"role": "user",   "content": json.dumps(names)},
            ],
            temperature=0,
            max_tokens=LLM_BATCH_SIZE * 30,
        )
        text = response.choices[0].message.content.strip()
        arr  = json.loads(text)
        if isinstance(arr, list) and len(arr) == len(names):
            return arr
    except Exception as e:
        print(f"  [groq batch error]: {e}")
    return None


def _call_ollama_batch(names: list[str], model: str) -> list[dict] | None:
    try:
        import requests
        response = requests.post(
            "http://localhost:11434/api/chat",
            json={
                "model":   model,
                "messages": [
                    {"role": "system", "content": LLM_BATCH_SYSTEM_PROMPT},
                    {"role": "user",   "content": json.dumps(names)},
                ],
                "stream":  False,
                "options": {"temperature": 0},
            },
            timeout=120,
        )
        response.raise_for_status()
        text = response.json()["message"]["content"].strip()
        text = re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL).strip()
        match = re.search(r"\[.*\]", text, re.DOTALL)
        if match:
            arr = json.loads(match.group())
            if isinstance(arr, list) and len(arr) == len(names):
                return arr
    except Exception as e:
        print(f"  [ollama batch error]: {e}")
    return None


def llm_normalise_batch(
    names: list[str],
    backend: str,
    model: str,
    cache: dict,
) -> dict[str, dict]:
    """
    Normalise a batch of gene names via LLM.
    Checks cache first, only calls LLM for uncached names.
    Returns { gene_name: {is_gene, normalised} } for all input names.
    """
    results   = {}
    to_lookup = []

    for name in names:
        if name in cache:
            results[name] = cache[name]
        else:
            to_lookup.append(name)

    if not to_lookup:
        return results

    # Process in sub-batches
    for i in range(0, len(to_lookup), LLM_BATCH_SIZE):
        batch = to_lookup[i : i + LLM_BATCH_SIZE]
        print(f"  [llm] batch {i // LLM_BATCH_SIZE + 1}: {len(batch)} genes via {backend}")

        if backend == "groq":
            arr = _call_groq_batch(batch, model)
        else:
            arr = _call_ollama_batch(batch, model)

        if arr:
            for name, item in zip(batch, arr):
                cache[name]   = item
                results[name] = item
        else:
            # On failure mark as unconverted so we don't lose them
            for name in batch:
                fallback      = {"is_gene": True, "normalised": name}
                cache[name]   = fallback
                results[name] = fallback

        time.sleep(LLM_RETRY_DELAY)

    save_llm_cache(cache)
    return results


# ── per-gene conversion ───────────────────────────────────────────────────────

def convert_gene(
    gene_obj:      dict,
    lookups:       dict,
    v21_mapping:   dict,
    funplantgenes: dict,
    llm_result:    dict | None = None,   # pre-computed {is_gene, normalised} or None
    llm_backend:   str  | None = None,   # for source label only
) -> dict:
    """
    Convert a single gene object.

    Stage 1 — direct chip lookup (TraesCS direct, funplantgenes, BioMart)
    Stage 2 — apply pre-computed LLM normalisation + retry chip lookup
    Stage 3 — not found → unconverted_gene (queued for BLAST)
    """
    gene_name = str(gene_obj.get("Gene", "")).strip()

    result = {
        "gene_input":        gene_name,
        "traescs_id":        None,
        "conversion_source": None,
        "tag":               None,
        "llm_normalised":    None,
        "pmid":              str(gene_obj.get("PMID", "")).strip(),
        "journal":           str(gene_obj.get("Journal", "")).strip(),
        "source_reference":  str(gene_obj.get("Source Reference", "")).strip(),
    }

    if not gene_name:
        result["tag"] = TAG_NOT_GENE
        return result

    # Stage 1 — direct chip lookup
    traescs, source = chip_lookup(gene_name, lookups, v21_mapping, funplantgenes)
    if traescs:
        result.update({
            "traescs_id":        traescs,
            "conversion_source": source,
            "tag":               TAG_CONVERTED,
        })
        return result

    # Stage 2 — apply LLM result if provided
    if llm_result is not None:
        if not llm_result.get("is_gene", True):
            result["tag"]            = TAG_NOT_GENE
            result["llm_normalised"] = "__not_a_gene__"
            return result

        normalised = (llm_result.get("normalised") or "").strip()
        if normalised and normalised != gene_name:
            result["llm_normalised"] = normalised
            traescs, source = chip_lookup(normalised, lookups, v21_mapping, funplantgenes)
            if traescs:
                result.update({
                    "traescs_id":        traescs,
                    "conversion_source": f"llm_{llm_backend or 'llm'}+{source}",
                    "tag":               TAG_CONVERTED,
                })
                return result

    # Stage 3 — queue for BLAST
    result["tag"] = TAG_UNCONVERTED_GENE
    return result


# ── main pipeline ─────────────────────────────────────────────────────────────

def run(llm_name: str, biomart_file: str, llm_backend: str | None, llm_model: str | None):
    gen_dir = Path(DEFAULT_GEN_BASE) / llm_name
    if not gen_dir.exists():
        print(f"Generation directory not found: {gen_dir}")
        return

    json_files = sorted(gen_dir.glob("*.json"))
    if Path(MAKER_FILE).exists():
        with open(MAKER_FILE, "r") as f:
            processed_genesets = {line.strip() for line in f if line.strip()}
        json_files = [f for f in json_files if f.stem in processed_genesets]
    print(f"\nFound {len(json_files)} gene set files in {gen_dir}\n")

    print("Building lookup tables...")
    lookups = build_lookup_tables(biomart_file)

    print("Loading v11→v21 mapping...")
    v21_mapping = load_v21_mapping()

    print("Loading funplantgenes mapping...")
    funplantgenes = load_funplantgenes_mapping()
    print()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    llm_cache: dict = {}
    if llm_backend and llm_model:
        llm_cache = load_llm_cache()
        print(f"LLM cache loaded: {len(llm_cache)} cached entries")

    all_results = []
    stats = {
        TAG_CONVERTED:        0,
        TAG_UNCONVERTED_GENE: 0,
        TAG_NOT_GENE:         0,
        "funplantgenes":      0,
        "llm_normalised":     0,
        "total_genes":        0,
        "total_gene_sets":    0
    }

    for i, json_file in enumerate(json_files, 1):
        gene_set_name = json_file.stem
        print(f"[{i}/{len(json_files)}] {gene_set_name}")

        try:
            with open(json_file, "r") as f:
                genes = json.load(f)
            if not isinstance(genes, list):
                continue
        except Exception as e:
            print(f"  Error reading {json_file.name}: {e}")
            continue

        stats["total_gene_sets"] += 1
        geneset_results = []

        # Collect genes needing LLM normalisation for this gene set
        # Eligibility: failed stage 1, has a PMID, not an obvious non-gene
        llm_candidates: dict[str, list[dict]] = {}   # gene_name → [gene_obj, ...]
        stage1_results: dict[str, dict] = {}

        for gene_obj in genes:
            gene_name = str(gene_obj.get("Gene", "")).strip()
            if not gene_name:
                continue
            traescs, source = chip_lookup(gene_name, lookups, v21_mapping, funplantgenes)
            if traescs:
                stage1_results[gene_name] = {"traescs": traescs, "source": source}
            elif llm_backend and llm_model:
                pmid = str(gene_obj.get("PMID", "")).strip()
                if pmid and not is_obvious_non_gene(gene_name):
                    llm_candidates.setdefault(gene_name, []).append(gene_obj)

        # Batch LLM call for all candidates in this gene set
        llm_results: dict[str, dict] = {}
        if llm_candidates and llm_backend and llm_model:
            unique_names = list(llm_candidates.keys())
            llm_results  = llm_normalise_batch(unique_names, llm_backend, llm_model, llm_cache)

        # Now convert each gene using pre-computed LLM results
        for gene_obj in genes:
            gene_name  = str(gene_obj.get("Gene", "")).strip()
            llm_result = llm_results.get(gene_name) if llm_backend else None
            result     = convert_gene(gene_obj, lookups, v21_mapping, funplantgenes, llm_result, llm_backend)
            result["gene_set"] = gene_set_name
            geneset_results.append(result)
            stats["total_genes"] += 1
            tag = result.get("tag", TAG_NOT_GENE)
            if tag in stats:
                stats[tag] += 1
            if result.get("conversion_source") == "funplantgenes":
                stats["funplantgenes"] += 1
            if result.get("llm_normalised") and result.get("tag") == TAG_CONVERTED:
                stats["llm_normalised"] += 1

        all_results.extend(geneset_results)

        # save incrementally every 50 gene sets
        if i % 50 == 0:
            partial = [
                r for r in all_results
                if r.get("gene_input", "").strip() and r.get("tag") != TAG_NOT_GENE
            ]
            with open(RESULTS_JSON, "w") as f:
                json.dump(partial, f, indent=2)

    # final save — strip not_gene and empty gene_input entries
    filtered_results = [
        r for r in all_results
        if r.get("gene_input", "").strip()
        and r.get("tag") != TAG_NOT_GENE
    ]
    with open(RESULTS_JSON, "w") as f:
        json.dump(filtered_results, f, indent=2)

    generate_viewer(RESULTS_JSON, VIEWER_HTML)

    print(f"\n{'─' * 50}")
    print(f"SUMMARY")
    print(f"  Gene sets processed   : {stats['total_gene_sets']}")
    print(f"  Total genes           : {stats['total_genes']}")
    print(f"  Converted (direct)    : {stats[TAG_CONVERTED]}")
    print(f"    of which funplantgenes : {stats['funplantgenes']}")
    print(f"    of which llm normalised: {stats['llm_normalised']}")
    print(f"  Unconverted (gene)    : {stats[TAG_UNCONVERTED_GENE]} ← queued for BLAST")
    print(f"  Not a gene            : {stats[TAG_NOT_GENE]}")
    print(f"\n  Results → {RESULTS_JSON}")
    print(f"  Viewer  → {VIEWER_HTML}")


# ── HTML viewer ───────────────────────────────────────────────────────────────

def generate_viewer(results_json: str, output_html: str):
    """
    Generate a lightweight paginated HTML viewer.
    Loads JSON data, renders 50 rows per page, supports search + tag filter.
    Only shows 'converted' and 'gene ⚑' flags — not_gene and empty entries
    are already stripped from the JSON before this is called.
    """
    html = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Gene conversion results</title>
<style>
  :root {
    --bg: #0f1117; --surface: #1a1d27; --border: #2a2d3a;
    --text: #e2e4ed; --muted: #8b8fa8; --accent: #6c8fff;
    --green: #4ade80; --amber: #fbbf24; --red: #f87171; --gray: #6b7280;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { background: var(--bg); color: var(--text); font-family: 'SF Mono', 'Fira Code', monospace; font-size: 13px; }
  header { padding: 20px 24px 16px; border-bottom: 1px solid var(--border); display: flex; align-items: center; gap: 16px; flex-wrap: wrap; }
  h1 { font-size: 15px; font-weight: 600; letter-spacing: .5px; color: var(--accent); flex: 1; }
  #stats { display: flex; gap: 12px; flex-wrap: wrap; }
  .stat { background: var(--surface); border: 1px solid var(--border); border-radius: 6px; padding: 4px 10px; font-size: 11px; }
  .stat span { font-weight: 600; }
  .controls { padding: 12px 24px; display: flex; gap: 10px; flex-wrap: wrap; border-bottom: 1px solid var(--border); }
  input[type=text] { background: var(--surface); border: 1px solid var(--border); color: var(--text); border-radius: 6px; padding: 6px 12px; font-size: 12px; font-family: inherit; width: 280px; }
  input[type=text]:focus { outline: none; border-color: var(--accent); }
  select { background: var(--surface); border: 1px solid var(--border); color: var(--text); border-radius: 6px; padding: 6px 10px; font-size: 12px; font-family: inherit; }
  .table-wrap { overflow-x: auto; }
  table { width: 100%; border-collapse: collapse; }
  th { background: var(--surface); color: var(--muted); font-size: 11px; font-weight: 500; text-align: left; padding: 10px 14px; border-bottom: 1px solid var(--border); white-space: nowrap; position: sticky; top: 0; }
  td { padding: 8px 14px; border-bottom: 1px solid var(--border); vertical-align: top; max-width: 300px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  tr:hover td { background: var(--surface); }
  .tag { display: inline-block; font-size: 10px; font-weight: 600; padding: 2px 7px; border-radius: 10px; white-space: nowrap; }
  .tag-converted            { background: #14532d; color: var(--green); }
  .tag-unconverted_gene     { background: #3d2a00; color: var(--amber); }
  .tag-gene                 { background: #3d2a00; color: var(--amber); }
  .traescs { font-size: 11px; color: var(--accent); }
  .pmid a { color: var(--muted); text-decoration: none; }
  .pmid a:hover { color: var(--accent); }
  .source { font-size: 10px; color: var(--muted); }
  .pagination { display: flex; align-items: center; gap: 8px; padding: 12px 24px; border-top: 1px solid var(--border); }
  .pagination button { background: var(--surface); border: 1px solid var(--border); color: var(--text); border-radius: 5px; padding: 5px 12px; cursor: pointer; font-family: inherit; font-size: 12px; }
  .pagination button:hover { border-color: var(--accent); color: var(--accent); }
  .pagination button:disabled { opacity: .3; cursor: default; }
  #page-info { color: var(--muted); font-size: 12px; margin: 0 8px; }
  .ref-cell { max-width: 250px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; cursor: help; }
  #loading { padding: 40px; text-align: center; color: var(--muted); }
</style>
</head>
<body>

<header>
  <h1>Gene conversion results</h1>
  <div id="stats"></div>
</header>

<div class="controls">
  <input type="text" id="search" placeholder="Search gene, gene set, TraesCS..." oninput="filterData()">
  <select id="tag-filter" onchange="filterData()">
    <option value="">All tags</option>
    <option value="converted">Converted ✓</option>
    <option value="unconverted_gene">Gene ⚑</option>
  </select>
  <select id="page-size" onchange="filterData()">
    <option value="50">50 per page</option>
    <option value="100">100 per page</option>
    <option value="200">200 per page</option>
  </select>
</div>

<div class="table-wrap">
  <div id="loading">Loading results...</div>
  <table id="table" style="display:none">
    <thead>
      <tr>
        <th>Tag</th>
        <th>Gene</th>
        <th>TraesCS ID (v2.1)</th>
        <th>Source</th>
        <th>Gene set</th>
        <th>PMID</th>
        <th>Journal</th>
        <th>Reference</th>
      </tr>
    </thead>
    <tbody id="tbody"></tbody>
  </table>
</div>

<div class="pagination">
  <button id="btn-prev" onclick="changePage(-1)" disabled>&#8592; Prev</button>
  <span id="page-info"></span>
  <button id="btn-next" onclick="changePage(1)">Next &#8594;</button>
  <span id="count-info" style="margin-left:auto;color:var(--muted);font-size:12px"></span>
</div>

<script>
let allData = [];
let filtered = [];
let page = 1;

const TAG_LABELS = {
  "converted":        { label: "converted ✓", cls: "tag-converted" },
  "unconverted_gene": { label: "gene ⚑",       cls: "tag-gene"      },
};

async function loadData() {
  try {
    const r = await fetch('conversion_results.json');
    allData = await r.json();
    renderStats();
    filterData();
    document.getElementById('loading').style.display = 'none';
    document.getElementById('table').style.display = 'table';
  } catch(e) {
    document.getElementById('loading').textContent = 'Error loading data: ' + e.message;
  }
}

function renderStats() {
  const counts = {};
  allData.forEach(r => { counts[r.tag] = (counts[r.tag] || 0) + 1; });
  const el = document.getElementById('stats');
  const labels = {
    converted:        ['Converted',       'green'],
    unconverted_gene: ['Queued for BLAST','amber'],
  };
  el.innerHTML = Object.entries(labels).map(([k, [label, color]]) =>
    `<div class="stat" style="color:var(--${color})"><span>${counts[k]||0}</span> ${label}</div>`
  ).join('') + `<div class="stat"><span>${allData.length}</span> total</div>`;
}

function filterData() {
  const q   = document.getElementById('search').value.toLowerCase();
  const tag = document.getElementById('tag-filter').value;

  filtered = allData.filter(r => {
    const matchTag = !tag || r.tag === tag;
    const matchQ   = !q || [r.gene_input, r.traescs_id,
                             r.gene_set, r.pmid, r.journal]
                            .some(v => v && v.toLowerCase().includes(q));
    return matchTag && matchQ;
  });

  page = 1;
  renderPage();
}

function renderPage() {
  const ps    = parseInt(document.getElementById('page-size').value);
  const start = (page - 1) * ps;
  const rows  = filtered.slice(start, start + ps);
  const tbody = document.getElementById('tbody');

  tbody.innerHTML = rows.map(r => {
    const tagInfo = TAG_LABELS[r.tag] || { label: r.tag, cls: 'tag-' + r.tag };
    return `
    <tr>
      <td><span class="tag ${tagInfo.cls}">${tagInfo.label}</span></td>
      <td>${esc(r.gene_input)}</td>
      <td>${r.traescs_id ? `<span class="traescs">${esc(r.traescs_id)}</span>` : '<span style="color:var(--muted)">—</span>'}</td>
      <td><span class="source">${esc(r.conversion_source || '—')}</span></td>
      <td style="max-width:200px;overflow:hidden;text-overflow:ellipsis">${esc(r.gene_set)}</td>
      <td class="pmid">${r.pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${r.pmid}" target="_blank">${r.pmid}</a>` : '—'}</td>
      <td>${esc(r.journal || '—')}</td>
      <td class="ref-cell" title="${esc(r.source_reference || '')}">${esc(r.source_reference || '—')}</td>
    </tr>`;
  }).join('');

  const total = filtered.length;
  const pages = Math.ceil(total / ps) || 1;
  document.getElementById('page-info').textContent = `Page ${page} of ${pages}`;
  document.getElementById('count-info').textContent = `${total.toLocaleString()} results`;
  document.getElementById('btn-prev').disabled = page <= 1;
  document.getElementById('btn-next').disabled = page >= pages;
}

function changePage(dir) {
  const ps    = parseInt(document.getElementById('page-size').value);
  const pages = Math.ceil(filtered.length / ps) || 1;
  page = Math.max(1, Math.min(pages, page + dir));
  renderPage();
  window.scrollTo(0, 0);
}

function esc(s) {
  if (!s) return '';
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}

loadData();
</script>
</body>
</html>"""

    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Viewer written → {output_html}")


# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Pre-BLAST gene conversion: BioMart + funplantgenes + optional LLM normalisation"
    )
    parser.add_argument(
        "--llm",
        required=True,
        help="LLM name used in generation e.g. deepseek-r1:8b"
    )
    parser.add_argument(
        "--biomart",
        default=DEFAULT_BIOMART,
        help=f"Merged BioMart TSV file (default: {DEFAULT_BIOMART})"
    )
    parser.add_argument(
        "--use-llm",
        action="store_true",
        help="Enable LLM normalisation for unconverted genes"
    )
    parser.add_argument(
        "--llm-backend",
        choices=["groq", "ollama"],
        default="ollama",
        help="LLM backend to use for normalisation (default: ollama)"
    )
    parser.add_argument(
        "--llm-model",
        default=None,
        help="Model name for LLM normalisation (default: llama3-8b-8192 for groq, deepseek-r1:8b for ollama)"
    )
    args = parser.parse_args()

    # Resolve default model per backend
    if args.use_llm and not args.llm_model:
        args.llm_model = DEFAULT_LLM_MODEL_GROQ if args.llm_backend == "groq" else DEFAULT_LLM_MODEL_OLLAMA

    llm_backend = args.llm_backend if args.use_llm else None
    llm_model   = args.llm_model   if args.use_llm else None

    if args.use_llm:
        print(f"LLM normalisation enabled: backend={llm_backend}, model={llm_model}")

    run(
        llm_name     = args.llm,
        biomart_file = args.biomart,
        llm_backend  = llm_backend,
        llm_model    = llm_model,
    )