import time
import json
import os
from langgraph.graph import END, StateGraph
from pubtator import Pubtator
from utils import GraphState, LLMGraphState, get_llm, get_llm_json_mode, clean_model_output, save_to_json_list
from langchain_core.messages import HumanMessage, SystemMessage
from instructs import rag_prompt2, grade_abstracts_instructions2, grade_full_texts_instructions

CHECKED_PMIDS_FILE = "checked_pmids_gene_checker.json"

def _build_relation_query(name: str, gene_ids: dict, include_geneset: bool = True) -> str:
    """
    Build a PubTator/grading relation query including both v2.1 and v1.1 IDs,
    transcript names, and gene names from both versions.
    include_geneset=True adds the gene set name + wheat context (for retrieval).
    include_geneset=False returns gene-only terms (for grading).
    """
    v21_id   = gene_ids.get("gene_stable_id", "")
    v11_id   = gene_ids.get("gene_stable_id_v11", "")
    tx_v11   = gene_ids.get("transcript_name", "")
    tx_v21   = gene_ids.get("transcript_name_v21", "")
    name_v11 = gene_ids.get("gene_name", "")
    name_v21 = gene_ids.get("gene_name_v21", "")

    gene_terms = [v21_id]

    if v11_id and v11_id != v21_id:
        gene_terms.append(v11_id)
    if tx_v21 and tx_v21 not in gene_terms:
        gene_terms.append(tx_v21)
    if tx_v11 and tx_v11 != tx_v21 and tx_v11 not in gene_terms:
        gene_terms.append(tx_v11)

    # Gene names — add whichever are present, non-nan, and not duplicates
    for gname in dict.fromkeys([name_v21, name_v11]):   # preserves order, deduplicates
        if gname and gname != "nan":
            gene_terms.append(f"@GENE_{gname}")

    gene_str = " OR ".join(filter(None, gene_terms))

    if include_geneset:
        return f"{name} in T aestivum OR Triticum aestivum OR wheat AND {gene_str}"
    return gene_str


def retrieve_pubtator_abstracts(state: LLMGraphState):
    """
    Retrieve abstracts for a given gene set.
    - If cached abstracts exist, load them from disk.
    - Otherwise, query PubTator, save the raw abstracts once, and return them.
    """
    geneset  = state["geneset"]
    name     = geneset["name"]
    gene_ids = geneset["gene_ids"]
    v21_id   = gene_ids["gene_stable_id"]
    print(f"Retrieving abstracts for gene set: {name} and gene: {v21_id}")

    base_dir   = "abstracts/gene_related_abstracts"
    os.makedirs(base_dir, exist_ok=True)
    cache_file = os.path.join(base_dir, f"{name}_{v21_id}.json")

    if os.path.exists(cache_file):
        print(f"Loading cached abstracts for {name} / {v21_id}...")
        with open(cache_file) as f:
            return {"documents": json.load(f)}

    relation_query = _build_relation_query(name, gene_ids, include_geneset=True)
    print(f"Querying PubTator with: {relation_query}")

    pmids     = Pubtator.search_pubtator_ID(relation=relation_query, limit=1)
    abstracts = []
    for pmid in pmids:
        try:
            abs_data = Pubtator.export_abstract(pmid, check_for_genes=False)
            if abs_data:
                abstracts.append(abs_data)
        except Exception as e:
            print(f"Error fetching PMID {pmid}: {e}")

    save_to_json_list(abstracts, cache_file)
    print(f"Saved {len(abstracts)} abstracts to {cache_file}")
    return {"documents": abstracts}

def retrieve_pubtator_full_text(state: LLMGraphState):
    """
    Retrieve full text for a given gene set.
    - If cached full text exist, load them from disk.
    - Otherwise, query PubTator, save the raw full text once, and return them.
    """
    geneset  = state["geneset"]
    name     = geneset["name"]
    gene_ids = geneset["gene_ids"]
    v21_id   = gene_ids["gene_stable_id"]
    print(f"Retrieving full text for gene set: {name} and gene: {v21_id}")

    base_dir   = "abstracts/gene_related_full_text"
    os.makedirs(base_dir, exist_ok=True)
    cache_file = os.path.join(base_dir, f"{name}_{v21_id}_full_text.json")

    if os.path.exists(cache_file):
        print(f"Loading cached full text for {name} / {v21_id}...")
        with open(cache_file) as f:
            return {"documents": json.load(f)}

    relation_query = _build_relation_query(name, gene_ids, include_geneset=True)
    print(f"Querying PubTator with: {relation_query}")

    pmids      = Pubtator.search_pubtator_ID(relation=relation_query, limit=1)
    full_texts = []
    for pmid in pmids:
        try:
            abs_data = Pubtator.export_full_text(pmid, check_for_genes=False)
            if abs_data:
                full_texts.append(abs_data)
        except Exception as e:
            print(f"Error fetching PMID {pmid}: {e}")

    save_to_json_list(full_texts, cache_file)
    print(f"Saved {len(full_texts)} full texts to {cache_file}")
    return {"documents": full_texts}


def grade_abstracts(state: LLMGraphState):
    """
    Grade abstracts for gene set relevance.
    Saves to state['documents_filtered'][llm_name].
    """
    llm_name = state.get("llm_name", "deepseek-r1:8b")
    geneset = state["geneset"]
    name = geneset["name"]
    gene_ids = geneset["gene_ids"]

    # 1. Access shared raw documents
    documents = state.get("documents", [])
    if not documents:
        print(f"No abstracts to grade for {name} and gene {gene_ids['gene_stable_id']}. Returning empty filtered list.")
        # Return empty list for this LLM
        return {"documents_filtered": {llm_name: []}}
    
    relation_query = _build_relation_query(name, gene_ids, include_geneset=False)
    print(f"Grading abstracts for {name} / {gene_ids['gene_stable_id']} using {llm_name}, query: {relation_query}...")

    question = (
        f"Does this abstract discuss BOTH the Triticum aestivum (common wheat) gene set "
        f"'{name}' meaning {geneset["definition"]} AND the gene {relation_query}? Include only wheat-relevant mechanisms and associations. "
        f"Answer 'yes' or 'no'."
    )

    llm = get_llm_json_mode(llm_name)
    filtered = []

    for doc in documents:
        abstract_text = (
            f"Title: {doc.get('title', '')}\n"
            f"Abstract: {doc.get('abstract', '')}"
        )
        
        try:
            result = llm.invoke([
                SystemMessage(content=grade_abstracts_instructions2),
                HumanMessage(content=f"Question: {question}\n\nAbstract:\n{abstract_text}")
            ])

            grade = json.loads(result.content)["binary_score"].strip().lower()
            if grade == "yes":
                filtered.append(doc)
        except Exception as e:
            print(f"[{llm_name}] Skipping abstract due to parse error: {e}")
            continue

    print(f"[{llm_name}] Kept {len(filtered)} abstracts after grading for {geneset}")
    
    # 2. Return using the nested structure
    return {"documents_filtered": {llm_name: filtered}}


def grade_full_texts(state: LLMGraphState):
    """
    Grade full texts for gene set relevance.
    Saves to state['documents_filtered'][llm_name].
    """
    llm_name = state.get("llm_name", "deepseek-r1:8b")
    geneset = state["geneset"]
    name = geneset["name"]
    gene_ids = geneset["gene_ids"]

    # 1. Access shared raw documents
    documents = state.get("documents", [])
    if not documents:
        print(f"No full texts to grade for {name} and gene {gene_ids['gene_stable_id']}. Returning empty filtered list.")
        # Return empty list for this LLM
        return {"documents_filtered": {llm_name: []}}
    
    relation_query = _build_relation_query(name, gene_ids, include_geneset=False)
    print(f"Grading full texts for {name} / {gene_ids['gene_stable_id']} using {llm_name}, query: {relation_query}...")

    question = (
        f"Does this full text section discuss BOTH the Triticum aestivum (common wheat) gene set "
        f"'{name}' meaning {geneset["definition"]} AND the gene {relation_query}? Include only wheat-relevant mechanisms and associations. "
        f"Answer 'yes' or 'no'."
    )

    llm = get_llm_json_mode(llm_name)
    filtered = []

    for doc in documents:
        # get only a section in the full text to grade, for example the results or discussion section, if available, otherwise fallback to the abstract
        doc_title = doc.get("title", "")
        doc_pmid = doc.get("pmid", "")
        sections = doc.get("sections", {})
        filtered_sections = {}

        for section_key, section in sections.items():
            section_text = (
                f"Section: {section.get('title', '')}\n"
                f"Text: {section.get('text', '')}"
            )

            try:
                result = llm.invoke([
                    SystemMessage(content=grade_abstracts_instructions2),
                    HumanMessage(content=f"Question: {question}\n\nFull Text Section:\n{section_text}")
                ])

                grade = json.loads(result.content)["binary_score"].strip().lower()
                if grade == "yes":
                    filtered_sections[section_key] = section
            except Exception as e:
                print(f"[{llm_name}] Skipping section {section_key} in full text for PMID {doc_pmid} due to parse error: {e}")
                continue
        
        if filtered_sections:
            filtered.append({
                "pmid": doc_pmid,
                "title": doc_title,
                "sections": filtered_sections
            })
        
    print(f"[{llm_name}] Kept {len(filtered)} full text sections after grading for {geneset}")
    
    # 2. Return using the nested structure
    return {"documents_filtered": {llm_name: filtered}}


def generate(state: LLMGraphState):
    """
    Use only the filtered abstracts for this specific LLM to validate the association.
    """
    llm_name = state.get("llm_name", "deepseek-r1:8b")

    geneset = state["geneset"]
    name = geneset["name"]
    gene_ids = geneset["gene_ids"]

    # 1. Retrieve filtered documents specifically for this LLM
    all_filtered = state.get("documents_filtered", {})
    documents = all_filtered.get(llm_name, [])

    if not documents:
        print(f"[{llm_name}] No filtered abstracts for {name} and gene {gene_ids['gene_stable_id']}")
        return {"generation": {llm_name: []}}

    pmids = [d.get("pmid") for d in documents]
    formatted_docs = [
        f"PMID: {d.get('pmid')}\nTitle: {d.get('title')}\nJournal: {d.get('journal')}\nAbstract: {d.get('abstract')}"
        for d in documents
    ]
    context = "\n\n".join(formatted_docs)

    relation_query = _build_relation_query(name, gene_ids, include_geneset=False)
    question = f"Is gene {relation_query} supported as being associated with Triticum aestivum (common wheat) gene set '{name}'?"

    llm = get_llm_json_mode(llm_name)
    
    try:
        result = llm.invoke([
            SystemMessage(content="You are a precise biomedical reasoning model specialized in Triticum aestivum (common wheat) gene set validation. Respond only in JSON."),
            HumanMessage(content=rag_prompt2.format(context=context, question=question))
        ])
        generation = json.loads(result.content)
    except Exception:
        try:
            generation = json.loads(clean_model_output(result.content))
        except:
            generation = {}

    generation["PMIDS"] = pmids

    # Save to disk — use v2.1 ID as filename
    outfile = f"out/geneset_checks/{llm_name}/{name}/{gene_ids['gene_stable_id']}.json"
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w") as f:
        json.dump(generation, f, indent=2)
    print(f"[{llm_name}] Saved generation for {geneset} to {outfile}")

    # 2. Return wrapped in a list to match LLMGraphState List type
    return {"generation": {llm_name: [generation]}}


def generate_full_texts(state: LLMGraphState):
    """
    Use only the filtered full text sections for this specific LLM to validate the association.
    """
    llm_name = state.get("llm_name", "deepseek-r1:8b")

    geneset = state["geneset"]
    name = geneset["name"]
    gene_ids = geneset["gene_ids"]

    # 1. Retrieve filtered documents specifically for this LLM
    all_filtered = state.get("documents_filtered", {})
    documents = all_filtered.get(llm_name, [])

    if not documents:
        print(f"[{llm_name}] No filtered abstracts for {name} and gene {gene_ids['gene_stable_id']}")
        return {"generation": {llm_name: []}}

    # run the generation for each full text document separately and save the results in a list, then aggregate the results at the end to return a single generation for the gene set, but also save the individual results for each full text section to disk for debugging and analysis
    pmids = [d.get("pmid") for d in documents]
    for doc in documents:
        doc_title = doc.get("title", "")
        doc_pmid = doc.get("pmid", "")
        sections = doc.get("sections", {})
        formatted_sections = [
            f"Section: {section.get('title', '')}\nText: {section.get('text', '')}"
            for section in sections.values()
        ]
        context = "\n\n".join(formatted_sections)

        relation_query = _build_relation_query(name, gene_ids, include_geneset=False)
        question = f"Is gene {relation_query} supported as being associated with Triticum aestivum (common wheat) gene set '{name}'?"

        llm = get_llm_json_mode(llm_name)
        try:
            result = llm.invoke([
                SystemMessage(content="You are a precise biomedical reasoning model specialized in Triticum aestivum (common wheat) gene set validation. Respond only in JSON."),
                HumanMessage(content=rag_prompt2.format(context=context, question=question))
            ])
            generation = json.loads(result.content)
        except Exception:
            # Fallback cleaning
            try:
                generation = json.loads(clean_model_output(result.content))
            except:
                generation = {}
        generation["PMIDS"] = [doc_pmid]

        # Save to disk
        outfile = f"out/geneset_checks/{llm_name}/{name}/{gene_ids['gene_stable_id']}.json"
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        with open(outfile, "w") as f:
            json.dump(generation, f, indent=2)
        print(f"[{llm_name}] Saved generation for {geneset} and PMID {doc_pmid} to {outfile}")

        # before checking another full text, check if the generation for the current full text supports the gene association, if yes, we can skip checking the remaining full texts to save time and resources, since we only need one piece of evidence to support the association, but if no, we continue checking the remaining full texts to see if there is any evidence in them
        if generation.get("Validation") == "yes":
            print(f"[{llm_name}] Found supporting evidence for {geneset} in PMID {doc_pmid}. Skipping remaining full texts.")
            break

    # 2. Return wrapped in a list to match LLMGraphState List type
    return {"generation": {llm_name: [generation]}}



def create_control_flow():
    # Ensure LLMGraphState in utils.py is updated to include `merge_dicts` logic!
    workflow = StateGraph(LLMGraphState)

    # 1. Retrieve Node
    workflow.add_node("retrieve", retrieve_pubtator_full_text)

    # 2. Grading Nodes (Parallel logic supported by StateGraph, though edges here are sequential)
    workflow.add_node("grade", lambda s: grade_full_texts(s))
    # workflow.add_node("grade_qwen", lambda s: grade_abstracts(s, "qwen3:32b"))
    # workflow.add_node("grade_deepseek", lambda s: grade_abstracts(s, "deepseek-r1:8b"))
    # workflow.add_node("grade_llama3", lambda s: grade_abstracts(s, "llama3.1:8b"))

    # 3. Generation Nodes
    workflow.add_node("generate", lambda s: generate_full_texts(s))
    # workflow.add_node("generate_qwen", lambda s: generate(s, "qwen3:32b"))
    # workflow.add_node("generate_deepseek", lambda s: generate(s, "deepseek-r1:8b"))
    # workflow.add_node("generate_llama3", lambda s: generate(s, "llama3.1:8b"))

    # Entry
    workflow.set_entry_point("retrieve")
    # # Flow
    # # Retrieve -> Grade Qwen -> Grade Deepseek -> Grade Llama
    # workflow.add_edge("retrieve", "grade_deepseek")
    # workflow.add_edge("grade_deepseek", "generate_deepseek")
    # workflow.add_edge("generate_deepseek",  END)


    # Flow
    # Retrieve -> Grade Qwen -> Grade Deepseek -> Grade Llama
    workflow.add_edge("retrieve", "grade")
    workflow.add_edge("grade", "generate")
    # workflow.add_edge("grade_qwen", "grade_deepseek")
    # workflow.add_edge("grade_deepseek", "grade_llama3")

    # # Grade Llama -> Gen Qwen -> Gen Deepseek -> Gen Llama
    # workflow.add_edge("grade_llama3", "generate_qwen")
    # workflow.add_edge("generate_qwen", "generate_deepseek")
    # workflow.add_edge("generate_deepseek", "generate_llama3")

    # End
    # workflow.add_edge("generate_llama3", END)
    workflow.add_edge("generate", END)

    return workflow.compile()

# def create_control_flow():
    # Ensure LLMGraphState in utils.py is updated to include `merge_dicts` logic!
    # workflow = StateGraph(LLMGraphState)

    # # 1. Retrieve Node
    # workflow.add_node("retrieve", retrieve_pubtator_abstracts)

    # # 2. Grading Nodes (Parallel logic supported by StateGraph, though edges here are sequential)
    # workflow.add_node("grade_qwen", lambda s: grade_abstracts(s, "qwen3:32b"))
    # workflow.add_node("grade_deepseek", lambda s: grade_abstracts(s, "deepseek-r1:8b"))
    # workflow.add_node("grade_llama3", lambda s: grade_abstracts(s, "llama3.1:8b"))

    # # 3. Generation Nodes
    # workflow.add_node("generate_qwen", lambda s: generate(s, "qwen3:32b"))
    # workflow.add_node("generate_deepseek", lambda s: generate(s, "deepseek-r1:8b"))
    # workflow.add_node("generate_llama3", lambda s: generate(s, "llama3.1:8b"))

    # # Entry
    # workflow.set_entry_point("retrieve")

    # # Flow
    # # Retrieve -> Grade Qwen -> Grade Deepseek -> Grade Llama
    # workflow.add_edge("retrieve", "grade_qwen")
    # workflow.add_edge("grade_qwen", "grade_deepseek")
    # workflow.add_edge("grade_deepseek", "grade_llama3")

    # # Grade Llama -> Gen Qwen -> Gen Deepseek -> Gen Llama
    # workflow.add_edge("grade_llama3", "generate_qwen")
    # workflow.add_edge("generate_qwen", "generate_deepseek")
    # workflow.add_edge("generate_deepseek", "generate_llama3")

    # # End
    # workflow.add_edge("generate_llama3", END)

    # return workflow.compile()