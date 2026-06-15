import os
import json
import argparse
import re

from utils import geneset_json_reader, phenotype_json_reader, read_gmt, read_geneset_to_gene_sets
from rag_pipeline_gene_set_maker_llm import create_control_flow as create_maker_flow
from rag_pipeline_gene_checker_llm import create_control_flow as create_checker_flow
from convert_gene_ids import build_id_lookup
from biomart_to_gmt import load_and_validate


GMT_PATH_GENE_NAME            = "out/wheat_gmt_gene_name.gmt"
GMT_PATH_STABLE_ID            = "out/wheat_gmt_gene_stable_id.gmt"
GMT_PATH_STABLE_ID_V21              = "out/wheat_gmt_gene_stable_id_2.1.gmt"
GMT_PATH_TRANSCRIPT_STABLE_ID       = "out/wheat_gmt_transcript_stable_id.gmt"
GMT_PATH_TRANSCRIPT_STABLE_ID_V21   = "out/wheat_gmt_transcript_stable_id_2.1.gmt"
BIOMART_EXPORT_TSV            = "data/mart_exports/mart_export.tsv"
BIOMART_EXPORT_TSV_V21        = "data/mart_exports/mart_export_2.1.tsv"
V11_TO_V21_MAPPING            = "v11_to_v21_mapping.json"

df     = load_and_validate(BIOMART_EXPORT_TSV)
lookup = build_id_lookup(df)

df_v21     = load_and_validate(BIOMART_EXPORT_TSV_V21)
lookup_v21 = build_id_lookup(df_v21)

# Build reverse mapping v2.1 → v1.1
with open(V11_TO_V21_MAPPING) as _f:
    _v11_to_v21 = json.load(_f)
v21_to_v11: dict[str, str] = {v: k for k, v in _v11_to_v21.items()}

def get_llm_file_paths(llm_name):
    """Generate LLM-specific file paths and ensure they exist."""
    os.makedirs("out", exist_ok=True)
    
    paths = {
        "processed_file": f"out/processed_genesets_{llm_name}.txt",
        "processed_maker_file": f"out/processed_genesets_maker_{llm_name}.txt",
        "processed_genes_file": f"out/processed_genes_{llm_name}.json",
        "processed_sets_file": f"out/processed_gene_sets_{llm_name}.txt"
    }
    
    # Create files if they don't exist
    for file_path in paths.values():
        if not os.path.exists(file_path):
            with open(file_path, "w") as f:
                if file_path.endswith(".json"):
                    f.write("{}")
                else:
                    f.write("")
    
    return paths


def load_processed(processed_file):
    """Load geneset names that have completed both pipelines."""
    if not os.path.exists(processed_file):
        return set()
    with open(processed_file, "r") as f:
        return set(line.strip() for line in f.readlines())


def mark_processed(geneset_name, processed_file):
    """Append a geneset to the processed file."""
    with open(processed_file, "a") as f:
        f.write(f"{geneset_name}\n")

def load_processed_genes(processed_genes_file):
    """Load processed genes per gene set for the checker pipeline."""
    if os.path.exists(processed_genes_file):
        with open(processed_genes_file, "r") as f:
            return json.load(f)
    return {}


def save_processed_genes(processed, processed_genes_file):
    """Save processed genes dictionary."""
    with open(processed_genes_file, "w") as f:
        json.dump(processed, f, indent=2)


def mark_gene_processed(gene_set, gene, processed, processed_genes_file):
    """Mark a gene as processed for a given gene set and persist to disk."""
    processed.setdefault(gene_set, [])
    if gene not in processed[gene_set]:
        processed[gene_set].append(gene)
        save_processed_genes(processed, processed_genes_file)


def load_completed_sets(processed_sets_file):
    """Load gene sets that have been completely checked."""
    if not os.path.exists(processed_sets_file):
        return set()
    with open(processed_sets_file, "r") as f:
        return set(line.strip() for line in f if line.strip())


def mark_set_complete(gene_set, processed_sets_file):
    """Append completed gene set to processed_gene_sets file."""
    with open(processed_sets_file, "a") as f:
        f.write(f"{gene_set}\n")


def run_checker_for_geneset(geneset, genes, llm_name, file_paths):
    """
    Run the checker pipeline for a single geneset name and its list of genes.
    Uses intersection logic: we only call this if the geneset exists in the GMT.
    genes are v2.1 transcript stable IDs from wheat_gmt_transcript_stable_id_2.1.gmt.
    """
    geneset_name = geneset["raw_go_name"].replace("/", " or ")
    print(f"Running checker pipeline for geneset: {geneset_name}")

    processed_genes = load_processed_genes(file_paths["processed_genes_file"])
    completed_sets  = load_completed_sets(file_paths["processed_sets_file"])

    if geneset_name in completed_sets:
        print(f"Gene set for {geneset_name} already completed. Skipping checker.")
        return

    processed_genes.setdefault(geneset_name, [])
    graph = create_checker_flow()

    for gene in genes:
        # gene is a v2.1 transcript stable ID — direct lookup works
        gene_ids = lookup_v21.get(gene)
        if not gene_ids:
            print(f"  Gene {gene} not found in v2.1 BioMart export.")
            continue

        gene_ids = dict(gene_ids)  # don't mutate the lookup

        # Attach v2.1 stable ID and find corresponding v1.1 ID
        v21_stable   = gene_ids.get("gene_stable_id", gene.split(".")[0])


        v11_stable   = v21_to_v11.get(re.sub(r"LC$", "", v21_stable, flags=re.I), "")
        gene_ids["gene_stable_id"]     = v21_stable
        gene_ids["gene_stable_id_v11"] = v11_stable

        if gene_ids["gene_stable_id"] in processed_genes.get(geneset_name, []):
            print(f"  Gene {gene} already processed — skipping")
            continue

        # Enrich with v1.1 gene name and transcript from v1.1 BioMart (for broader search)
        v11_info = lookup.get(v11_stable) if v11_stable else None
        gene_ids["transcript_name_v21"] = gene_ids.get("transcript_name", "")
        gene_ids["gene_name_v21"]       = gene_ids.get("gene_name", "")
        if v11_info:
            gene_ids["transcript_name"] = v11_info.get("transcript_name", "")
            gene_ids["gene_name"]       = v11_info.get("gene_name", "")

        geneset_state = {
            "name":       geneset_name,
            "gene_ids":   gene_ids,
            "definition": geneset["go_term_definition"]
        }


        print(f"  Checking {gene} for {geneset_name}")

        try:
            for _ in graph.stream({"geneset": geneset_state, "llm_name": llm_name}, stream_mode="values"):
                pass
            mark_gene_processed(geneset_name, gene, processed_genes, file_paths["processed_genes_file"])
            print(f"  Completed {gene} for {geneset_name}")
        except Exception as e:
            print(f"  Error processing {gene} in {geneset_name}: {e}")

    mark_set_complete(geneset_name, file_paths["processed_sets_file"])
    print(f"Completed checker for gene set: {geneset_name}")


def main():
    parser = argparse.ArgumentParser(description="Run maker and checker pipelines for genesets.")
    parser.add_argument(
        "--input_file",
        type=str,
        default="out/go_terms.json",
        help="Path to the input JSON file with geneset details."
    )
    # add llm to use as argument
    parser.add_argument(
        "--llm",
        type=str,
        default="deepseek-r1:8b",
        help="LLM to use for grading abstracts."
    )
    args = parser.parse_args()

    llm_name = args.llm
    # Generate LLM-specific file paths
    file_paths = get_llm_file_paths(llm_name)

    # Load genesets
    genesets = phenotype_json_reader(args.input_file)
    # print some of the loaded genesets for debugging
    print(f"Loaded {len(genesets)} genesets. Sample:")
    for gs in genesets[:3]:
        print(gs)

    # Load already fully processed genesets (maker + checker)
    processed = load_processed(file_paths["processed_file"])

    # Load GMT gene sets once
    if not os.path.exists(GMT_PATH_STABLE_ID_V21):
        raise FileNotFoundError(f"v2.1 GMT file not found at {GMT_PATH_STABLE_ID_V21}")
    gene_sets_gene_name                = read_geneset_to_gene_sets(GMT_PATH_GENE_NAME)
    gene_sets_stable_id                = read_geneset_to_gene_sets(GMT_PATH_STABLE_ID)
    gene_sets_stable_id_v21            = read_geneset_to_gene_sets(GMT_PATH_STABLE_ID_V21)
    gene_sets_transcript_stable_id     = read_geneset_to_gene_sets(GMT_PATH_TRANSCRIPT_STABLE_ID)
    gene_sets_transcript_stable_id_v21 = read_geneset_to_gene_sets(GMT_PATH_TRANSCRIPT_STABLE_ID_V21)
    # print all the geneset names in the GMT for debugging
    # print(f"Gene sets in GMT (transcript stable ID): {list(gene_sets_transcript_stable_id.keys())[:200]}")
    # exit(0)

    # Intersection logic:
    #   Checker only runs for genesets whose name appears in gene_sets.
    to_process = [p for p in genesets if p["raw_go_name"] not in processed]

    print(f"Total genesets in file: {len(genesets)}")
    print(f"Already fully processed (maker + checker): {len(processed)}")
    print(f"Remaining to process: {len(to_process)}")
    print(f"Using LLM-specific files with suffix: {llm_name}")

    if not to_process:
        print("All genesets already processed. Nothing to do.")
        return


    for geneset in to_process:
        maker_complete = False
        checker_complete = False
        name = geneset["raw_go_name"]
        print(f"\nProcessing geneset: {name}")
        # replace "/" with "or" in name to avoid file path issues
        name = name.replace("/", " or ")

        # Maker pipeline
        try:
            # check if this geneset has already completed the maker pipeline by looking for its name in processed_maker_file
            processed_maker = load_processed(file_paths["processed_maker_file"])
            if name in processed_maker:
                print(f"Maker pipeline already completed for {name}. Skipping maker.")
                maker_complete = True
            else:
                maker_graph = create_maker_flow()
                inputs = {"geneset": geneset, "llm_name": llm_name}

                for _ in maker_graph.stream(inputs, stream_mode="values"):
                    pass

                print(f"Maker pipeline completed for {name}")
                # mark this geneset as having completed the maker pipeline by adding to processed_maker_file
                mark_processed(name, file_paths["processed_maker_file"])
                maker_complete = True
        except Exception as e:
            print(f"Error in maker pipeline for {name}: {e}")
            # Do not mark as processed; continue to next geneset
            continue

        # Checker pipeline only if this geneset appears in one of the GMT files
        # make name and all gene_sets keys lowercase for case-insensitive comparison
        # name = geneset["definition"].lower()
        
        
        if name.capitalize().replace("-", " ") in gene_sets_transcript_stable_id_v21:
            genes = gene_sets_transcript_stable_id_v21[name.capitalize().replace("-", " ")]
            run_checker_for_geneset(geneset, genes, llm_name, file_paths)
            checker_complete = True
        else:
            print(f"No matching gene set in gene set database for geneset '{name}'. Skipping checker for this geneset.")

        # Only now mark geneset as fully processed
        if maker_complete and checker_complete:
            mark_processed(name, file_paths["processed_file"])
            print(f"Finished geneset: {name}")

    print(f"\nAll genesets processed. Progress saved.")


if __name__ == "__main__":
    main()