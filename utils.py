import operator
from typing_extensions import TypedDict
from typing import Annotated, Dict, List
### LLM
from langchain_ollama import ChatOllama
import re
import json
import mygene
import os
import csv
from collections import defaultdict

# legacy graph state
# class GraphState(TypedDict):
#     """
#     Graph state is a dictionary that contains information we want to propagate to, and modify in, each graph node.
#     """

#     question: str  # User question
#     generation: str  # LLM generation
#     max_retries: int  # Max number of retries for answer generation
#     answers: int  # Number of answers generated
#     loop_step: Annotated[int, operator.add]
#     documents: List[str]  # List of retrieved documents

# 1. Define the reducer function
def merge_dicts(left: Dict, right: Dict) -> Dict:
    """
    Merges two dictionaries. 
    If keys overlap, the 'right' (new) value overwrites the 'left' (old).
    Example: {A:1} + {B:2} -> {A:1, B:2}
    """
    if not left:
        left = {}
    if not right:
        right = {}
    return {**left, **right}

class GraphState(TypedDict):
    geneset: str 
    documents: list # Raw documents
    # This allows multiple LLM nodes to add their specific results 
    # to this dictionary without deleting others.
    documents_filtered: Annotated[Dict[str, List], merge_dicts]
    
    # You likely want the same pattern for generation results
    generation: Annotated[Dict[str, List], merge_dicts]

class LLMGraphState(TypedDict):
    geneset: str 
    documents: list # Raw documents
    documents_fulltext: list
    documents_filtered: list
    documents_fulltext_filtered: list
    llm_name: str
    generation: list

class AAGraphState(TypedDict):
    geneset: str 
    num_of_abstracts: int
    documents: list # Raw documents
    documents_filtered: list
    llm_name: str
    generation: list

class DPGraphState(TypedDict):
    geneset: str
    llm_name: str
    generation: list

# Post-processing
def format_docs(docs):
    return "\n\n".join(doc.page_content for doc in docs)

def format_json_docs(docs):
    pass
 

# local_llm = "llama3.1:8b"
# llm = ChatOllama(model=local_llm, temperature=0)
# llm_json_mode = ChatOllama(model=local_llm, temperature=0, format="json")

def get_llm(local_llm="llama3.1:8b"):
    llm = ChatOllama(model=local_llm, temperature=0)
    return llm

def get_llm_json_mode(local_llm="llama3.1:8b"):
    llm_json_mode = ChatOllama(model=local_llm, temperature=0, format="json")
    return llm_json_mode

def write_gmt(file_path, gene_sets):
    with open(file_path, "a") as file:
        for gene_set, genes in gene_sets.items():
            entrez_genes, _, _ = id_mapping(genes)
            entrez_genes = list(set(entrez_genes))
            line = f"{gene_set}\tNA\t" + "\t".join(entrez_genes) + "\n"
            file.write(line)

def load_gmt(filepath):
    """Load a GMT file into a dict: {gene_set_name: set(genes)}"""
    gene_sets = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 2:
                gene_sets[parts[0]] = set(parts[2:])
    return gene_sets

def id_mapping(genes, mode='entrezgene'):
    mg = mygene.MyGeneInfo()
    out = []
    if mode =='entrezgene':
        out = mg.querymany(genes, scopes='symbol,reporter,accession,entrezgene', fields='entrezgene', species='human')
    elif mode == 'symbol':
        out = mg.querymany(genes, scopes='symbol,reporter,accession,entrezgene', fields='symbol', species='human')
    valid_genes = []
    mapped_genes = []
    invalid_genes = []
    for gene_info in out:
        if "notfound" in gene_info:
            # print("got here")
            invalid_genes.append(gene_info["query"])
        else:
            valid_genes.append(gene_info["query"])
            if 'entrezgene' in gene_info:
                mapped_genes.append(gene_info[mode])
    return mapped_genes, valid_genes, invalid_genes

def parse_out_json(content):
    # Use regex to capture text between content= and additional_kwargs
    match = re.search(r"content=(.*?)additional_kwargs", content, re.DOTALL)
    
    if match:
        # Extract the JSON-like string
        json_str = match.group(1).strip()
        
        # Clean up the JSON string
        json_str = json_str.replace('\\n', '')  # Remove newline escape
        json_str = json_str.replace('\\r', '')  # Remove carriage returns
        json_str = json_str.replace('\\t', '')  # Remove tab escape
        # json_str = json_str.replace('\\"', '')  # Fix escaped double quotes
        json_str = json_str.replace("'", '')  # Fix escaped single quotes
        json_str = json_str.replace('\\\\', '')    # Replace double backslash with single
        json_str = json_str.replace('\\', '')   

        # Optional: Remove triple backticks (if any)
        if json_str.startswith('```') and json_str.endswith('```'):
            json_str = json_str[3:-3].strip()

        print("Cleaned JSON string:", json_str)

        # Try to parse the cleaned string as JSON
        try:
            parsed = json.loads(json_str)
            return parsed
        except json.JSONDecodeError as e:
            print("JSON decoding failed:", e)
            return None
    else:
        print("No match found for the pattern")
        return None
    
import json

def geneset_json_reader(json_file, attributes=["Standard Gene Set Name","Full Description/Abstract"]):
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)

    result = []
    if isinstance(data, list):
        for item in data:
            # get the "Standard Gene Set Name" as "name" and "Full Description/Abstract" as "description"
            new_entry = {
                "name": item.get("Standard Gene Set Name", None),
                "definition": " ".join(item.get("Full Description/Abstract", "").split(" ")[1:]) if item.get("Full Description/Abstract", None) else None
            }
            result.append(new_entry)
    elif isinstance(data, dict):
        print(f"Loaded dict with {len(data)} keys from {json_file}")
        for key, item in data.items():
            new_entry = {
                "name": item.get("Standard Gene Set Name", None),
                "definition": " ".join(item.get("Full Description/Abstract", "").split(" ")[1:]) if item.get("Full Description/Abstract", None) else None
            }
            result.append(new_entry)
    else:
        raise ValueError("JSON format not supported. Must be dict or list of dicts.")

    return result



def phenotype_json_reader(json_file, attributes=["go_id","raw_go_name", "go_term_definition", "go_domain"]):
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)

    result = []
    if isinstance(data, list):
        for item in data:
            result.append({attr: item.get(attr, None) for attr in attributes})
    elif isinstance(data, dict):
        result.append({attr: data.get(attr, None) for attr in attributes})
    else:
        raise ValueError("JSON format not supported. Must be dict or list of dicts.")

    return result


def clean_model_output(raw_output: str):
    # Remove <think>...</think> sections
    cleaned = re.sub(r"<think>.*?</think>", "", raw_output, flags=re.DOTALL)
    # Also remove any stray markdown or code fences
    cleaned = cleaned.strip().replace("```json", "").replace("```", "").strip()
    return cleaned




# Reads a 'phenotypes-to-genes' txt file and creates a phenotype-to-gene-sets file.
def build_phenotype_to_gene_sets(input_file: str, output_file: str):
    
    phenotype_genes = defaultdict(set)
    phenotype_names = {}

    # --- Read input file ---
    with open(input_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            hpo_id = row["hpo_id"].strip()
            hpo_name = row["hpo_name"].strip()
            gene_symbol = row["gene_symbol"].strip()

            phenotype_names[hpo_id] = hpo_name
            phenotype_genes[hpo_id].add(gene_symbol)

    # --- Write output file ---
    with open(output_file, "w", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["hpo_id", "hpo_name", "genes"])

        for hpo_id, genes in phenotype_genes.items():
            writer.writerow([hpo_id, phenotype_names[hpo_id], ",".join(sorted(genes))])


def compare_to_phenotypes_msigdb(phenotype_file, hpo_db_file):
    # --- Load phenotype gene sets ---
    phenotype_names = set()
    phenotype_ids = {}

    with open(phenotype_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            hpo_id = row["hpo_id"].strip()
            hpo_name = row["hpo_name"].strip()
            hpo_name = "HP_" + hpo_name
            hpo_name = hpo_name.replace("-", "_")
            phenotype_names.add(hpo_name.upper().replace(" ", "_"))  # normalize
            phenotype_ids[hpo_id] = hpo_name



    # --- Load MSigDB HPO gene set database ---
    db_gene_sets = set()
    db_full = {}

    with open(hpo_db_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if not parts:
                continue
            gs_name = parts[0].strip()
            # gs_name = parts[0].strip('HP_')
            db_gene_sets.add(gs_name)
            db_full[gs_name] = parts[1:]  # keep the rest if needed
    print(list(phenotype_names)[:2])

    print(list(db_gene_sets)[:2])

    # --- Compare ---
    in_both = phenotype_names & db_gene_sets
    only_in_phenotypes = phenotype_names - db_gene_sets
    only_in_db = db_gene_sets - phenotype_names

    print(f"Matches found: {len(in_both)}")
    print(f"Only in phenotypes file: {len(only_in_phenotypes)}")
    print(f"Only in database file: {len(only_in_db)}")

    return in_both, only_in_phenotypes, only_in_db, db_gene_sets

def check_is_gene_annotated(pmids, ga_pmids):
    # if not os.path.exists(PMIDS_FILE):
    #     raise FileNotFoundError(f"PMID list file {PMIDS_FILE} not found.")

    # Use a set for fast lookup
  

    # print(type(pmids[1]))
    # print(type(ga_pmids[1]))
    # Filter input pmids
    return [pid for pid in pmids if pid in ga_pmids]

    
def save_to_json_list(results, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)


def read_gmt(file_path):
    """
    Read a GMT file and return a dictionary:
    { gene_set_name: [list_of_genes] }

    Cleans gene set names such as 'HP_11_PAIRS_OF_RIBS' → '11 pairs of ribs'.
    This must match the phenotype['name'] used in the maker pipeline.
    """
    gene_sets = {}
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            raw_name = parts[0]
            cleaned_name = (
                raw_name.replace("HP_", "")
                        .replace("MP_", "")
                        .replace("_", " ")
                        .strip()
                        .capitalize()
            )
            genes = parts[2:]
            gene_sets[cleaned_name] = genes
    return gene_sets

# read plant gmt
def read_geneset_to_gene_sets(file_path):
    gene_sets = {}
    # there is no header in the plant gmt file, so we can read directly
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            raw_name = parts[0]
            cleaned_name = (
                raw_name.replace("HP_", "")
                        .replace("MP_", "")
                        .replace("_", " ")
                        .strip()
                        .capitalize()
            )
            genes = parts[2:]
            gene_sets[cleaned_name] = genes
    return gene_sets