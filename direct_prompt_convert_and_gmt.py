"""
direct_prompt_convert_and_gmt.py
---------------------------------
Runs gene conversion on the direct-prompting generation folder, then
curates a GMT file from the converted results.

Direct-prompting generations live at:
    out/geneset_generations_direct/<model>/

Each file is a JSON array of gene objects:
    [{"Gene": "...", "Source Reference": "...", "PMID": "", "Journal": ""}, ...]

Steps:
  1. Read all JSON files from the direct-prompting folder for the given model.
  2. Run the same BioMart + funplantgenes + optional LLM conversion logic
     from gene_convert.py.
  3. Write conversion_results to out/conversion_direct/<model>/conversion_results.json
     and a viewer HTML alongside it.
  4. Build a GMT file at out/conversion_direct/<model>/direct_prompt_curated.gmt
     containing only converted genes, one row per gene set.

Usage:
    python direct_prompt_convert_and_gmt.py --model deepseek-r1:8b
    python direct_prompt_convert_and_gmt.py --model deepseek-r1:8b --use-llm --llm-backend ollama
    python direct_prompt_convert_and_gmt.py --model deepseek-r1:8b --biomart out/merged_biomart_with_v11.tsv
"""

import argparse
import json
import os
import sys
import re
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Import shared logic from gene_convert.py (must be in same directory)
# ---------------------------------------------------------------------------
try:
    import gene_convert as gc
except ImportError:
    print("ERROR: gene_convert.py not found in the current directory.", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DIRECT_GEN_BASE  = "out/geneset_generations_direct"
DIRECT_OUT_BASE  = "out/conversion_direct"


# ---------------------------------------------------------------------------
# GMT builder
# ---------------------------------------------------------------------------

def build_gmt(conversion_results: list[dict], out_path: Path, model: str) -> None:
    """
    Build a GMT file from converted genes only.

    GMT format (tab-separated):
        <gene_set_name>  <description>  <gene1>  <gene2>  ...

    Description column: "direct_prompt:<model>"
    Only entries with tag == 'converted' and a non-empty traescs_id are included.
    Gene sets with zero converted genes are omitted.
    """
    from collections import defaultdict

    gene_sets: dict[str, list[str]] = defaultdict(list)
    seen: dict[str, set] = defaultdict(set)  # deduplicate per gene set

    for r in conversion_results:
        if r.get("tag") != gc.TAG_CONVERTED:
            continue
        traescs = (r.get("traescs_id") or "").strip()
        gene_set = (r.get("gene_set") or "").strip()
        if not traescs or not gene_set:
            continue
        if traescs not in seen[gene_set]:
            gene_sets[gene_set].append(traescs)
            seen[gene_set].add(traescs)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    description = f"direct_prompt:{model}"
    written = 0

    with out_path.open("w", encoding="utf-8") as f:
        for gene_set in sorted(gene_sets):
            genes = gene_sets[gene_set]
            if not genes:
                continue
            cols = [gene_set, description] + genes
            f.write("\t".join(cols) + "\n")
            written += 1

    print(f"\nGMT written → {out_path}")
    print(f"  Gene sets with converted genes : {written}")
    total_genes = sum(len(v) for v in gene_sets.values())
    print(f"  Total unique TraesCS IDs       : {total_genes}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_direct(
    model: str,
    biomart_file: str,
    llm_backend: str | None,
    llm_model: str | None,
) -> None:
    gen_dir = Path(DIRECT_GEN_BASE) / model
    if not gen_dir.exists():
        print(f"ERROR: Generation directory not found: {gen_dir}", file=sys.stderr)
        sys.exit(1)

    json_files = sorted(gen_dir.glob("*.json"))
    print(f"Found {len(json_files)} gene set files in {gen_dir}\n")

    if not json_files:
        print("Nothing to process.")
        return

    # Output paths for this model
    out_dir      = Path(DIRECT_OUT_BASE) / model
    results_json = out_dir / "conversion_results.json"
    viewer_html  = out_dir / "viewer.html"
    gmt_path     = out_dir / "direct_prompt_curated.gmt"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Override cache path so it doesn't collide with main pipeline cache
    gc.LLM_CACHE_FILE = str(out_dir / "llm_normalise_cache.json")

    # Build lookup tables
    print("Building lookup tables...")
    lookups = gc.build_lookup_tables(biomart_file)

    print("Loading v11→v21 mapping...")
    v21_mapping = gc.load_v21_mapping()

    print("Loading funplantgenes mapping...")
    funplantgenes = gc.load_funplantgenes_mapping()
    print()

    llm_cache: dict = {}
    if llm_backend and llm_model:
        llm_cache = gc.load_llm_cache()
        print(f"LLM cache loaded: {len(llm_cache)} cached entries")

    all_results = []
    stats = {
        gc.TAG_CONVERTED:        0,
        gc.TAG_UNCONVERTED_GENE: 0,
        gc.TAG_NOT_GENE:         0,
        "funplantgenes":         0,
        "llm_normalised":        0,
        "total_genes":           0,
        "total_gene_sets":       0,
    }

    for i, json_file in enumerate(json_files, 1):
        gene_set_name = json_file.stem
        print(f"[{i}/{len(json_files)}] {gene_set_name}")

        try:
            with open(json_file, "r", encoding="utf-8") as f:
                genes = json.load(f)
            if not isinstance(genes, list):
                print(f"  Skipping — not a list")
                continue
        except Exception as e:
            print(f"  Error reading {json_file.name}: {e}")
            continue

        stats["total_gene_sets"] += 1

        # Stage 1 for all genes + collect LLM candidates
        llm_candidates: dict[str, list[dict]] = {}
        stage1_hits: dict[str, dict] = {}

        for gene_obj in genes:
            gene_name = str(gene_obj.get("Gene", "")).strip()
            if not gene_name:
                continue
            traescs, source = gc.chip_lookup(gene_name, lookups, v21_mapping, funplantgenes)
            if traescs:
                stage1_hits[gene_name] = {"traescs": traescs, "source": source}
            elif llm_backend and llm_model:
                # Direct-prompt files have no PMIDs, so skip PMID gate
                if not gc.is_obvious_non_gene(gene_name):
                    llm_candidates.setdefault(gene_name, []).append(gene_obj)

        # Batch LLM normalisation
        llm_results: dict[str, dict] = {}
        if llm_candidates and llm_backend and llm_model:
            unique_names = list(llm_candidates.keys())
            llm_results  = gc.llm_normalise_batch(unique_names, llm_backend, llm_model, llm_cache)

        # Convert each gene
        for gene_obj in genes:
            gene_name  = str(gene_obj.get("Gene", "")).strip()
            llm_result = llm_results.get(gene_name) if llm_backend else None
            result     = gc.convert_gene(gene_obj, lookups, v21_mapping, funplantgenes, llm_result, llm_backend)
            result["gene_set"] = gene_set_name
            all_results.append(result)
            stats["total_genes"] += 1
            tag = result.get("tag", gc.TAG_NOT_GENE)
            if tag in stats:
                stats[tag] += 1
            if result.get("conversion_source") == "funplantgenes":
                stats["funplantgenes"] += 1
            if result.get("llm_normalised") and result.get("tag") == gc.TAG_CONVERTED:
                stats["llm_normalised"] += 1

        # Incremental save every 50 gene sets
        if i % 50 == 0:
            partial = [
                r for r in all_results
                if r.get("gene_input", "").strip() and r.get("tag") != gc.TAG_NOT_GENE
            ]
            with open(results_json, "w", encoding="utf-8") as f:
                json.dump(partial, f, indent=2)

    # Final save
    filtered_results = [
        r for r in all_results
        if r.get("gene_input", "").strip() and r.get("tag") != gc.TAG_NOT_GENE
    ]
    with open(results_json, "w", encoding="utf-8") as f:
        json.dump(filtered_results, f, indent=2)
    print(f"\nConversion results → {results_json}")

    # Viewer
    gc.generate_viewer(str(results_json), str(viewer_html))

    # GMT
    build_gmt(filtered_results, gmt_path, model)

    # Summary
    print(f"\n{'─' * 50}")
    print(f"SUMMARY — direct prompt: {model}")
    print(f"  Gene sets processed        : {stats['total_gene_sets']}")
    print(f"  Total genes                : {stats['total_genes']}")
    print(f"  Converted                  : {stats[gc.TAG_CONVERTED]}")
    print(f"    of which funplantgenes   : {stats['funplantgenes']}")
    print(f"    of which llm normalised  : {stats['llm_normalised']}")
    print(f"  Unconverted (gene)         : {stats[gc.TAG_UNCONVERTED_GENE]}")
    print(f"  Not a gene                 : {stats[gc.TAG_NOT_GENE]}")
    print(f"\n  Results  → {results_json}")
    print(f"  Viewer   → {viewer_html}")
    print(f"  GMT      → {gmt_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert direct-prompt gene set JSONs and build a curated GMT."
    )
    parser.add_argument(
        "--model",
        required=True,
        help="Model subfolder name e.g. deepseek-r1:8b"
    )
    parser.add_argument(
        "--biomart",
        default=gc.DEFAULT_BIOMART,
        help=f"Merged BioMart TSV (default: {gc.DEFAULT_BIOMART})"
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
        help="LLM backend (default: ollama)"
    )
    parser.add_argument(
        "--llm-model",
        default=None,
        help="Model for LLM normalisation (default: deepseek-r1:8b for ollama, llama3-8b-8192 for groq)"
    )
    args = parser.parse_args()

    if args.use_llm and not args.llm_model:
        args.llm_model = (
            gc.DEFAULT_LLM_MODEL_GROQ if args.llm_backend == "groq"
            else gc.DEFAULT_LLM_MODEL_OLLAMA
        )

    llm_backend = args.llm_backend if args.use_llm else None
    llm_model   = args.llm_model   if args.use_llm else None

    if args.use_llm:
        print(f"LLM normalisation: backend={llm_backend}, model={llm_model}")

    run_direct(
        model        = args.model,
        biomart_file = args.biomart,
        llm_backend  = llm_backend,
        llm_model    = llm_model,
    )
