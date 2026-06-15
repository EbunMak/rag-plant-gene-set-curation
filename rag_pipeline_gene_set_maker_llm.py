import time
import json
import os
import re
from langgraph.graph import END, StateGraph
from pubtator import Pubtator
from utils import LLMGraphState, get_llm, get_llm_json_mode, clean_model_output, check_is_gene_annotated
from langchain_core.messages import HumanMessage, SystemMessage
from instructs import rag_prompt, grade_abstracts_instructions, grade_full_texts_instructions, build_grader_question, abstract_grader_prompt

CHECKED_PMIDS_FILE = "checked_pmids.json"
PMIDS_FILE = "abstracts/pmids.txt"
VERIFY_MODEL = "qwen3:32b"   # model used for metadata verification
FULL_TEXT_LIMIT = 50         # max full texts to fetch per gene set

# Load gene-annotated PMIDs once
if os.path.exists(PMIDS_FILE):
    with open(PMIDS_FILE, "r") as f:
        ga_pmids = list({int(line.strip()) for line in f if line.strip()})
else:
    ga_pmids = []
    print(f"Warning: {PMIDS_FILE} not found. No PMIDs loaded.")


# ── helpers ───────────────────────────────────────────────────────────────────

def append_safe_save(new_genes: list, json_outfile: str) -> int:
    """
    Append new gene objects to an existing JSON file, skipping duplicates.
    Deduplication key: (Gene, PMID) tuple — full object must be unique on both.
    Returns the number of new genes actually written.
    """
    existing = []
    if os.path.exists(json_outfile):
        try:
            with open(json_outfile, "r") as f:
                existing = json.load(f)
            if not isinstance(existing, list):
                existing = []
        except Exception:
            existing = []

    # Build dedup set from existing entries
    seen = {
        (str(g.get("Gene", "")).strip().lower(), str(g.get("PMID", "")).strip())
        for g in existing
    }

    added = 0
    for gene_obj in new_genes:
        print(f"Processing gene object: {gene_obj}")
        key = (
            str(gene_obj.get("Gene", "")).strip().lower(),
            str(gene_obj.get("PMID", "")).strip()
        )
        if key not in seen:
            existing.append(gene_obj)
            seen.add(key)
            added += 1

    with open(json_outfile, "w") as f:
        json.dump(existing, f, indent=2)

    return added


def repair_json_with_regex(raw: str) -> str:
    """
    Stage 1: cheap regex-based cleanup before attempting an LLM repair call.
    Handles the most common failure modes across DeepSeek R1, Qwen, and Llama.
    """
    # Fix 1 — strip DeepSeek R1 <think>...</think> reasoning blocks
    raw = re.sub(r"<think>.*?</think>", "", raw, flags=re.DOTALL).strip()

    # Fix 2 — extract JSON from markdown code fences ```json ... ``` or ``` ... ```
    match = re.search(r"```(?:json)?\s*([\s\S]*?)```", raw)
    if match:
        raw = match.group(1).strip()

    # Fix 3 — find the first [ ... ] array anywhere in the output
    match = re.search(r"(\[[\s\S]*\])", raw)
    if match:
        raw = match.group(1).strip()

    return raw


def repair_json_with_llm(raw: str, repair_llm_name: str = "llama3.1:8b") -> list:
    """
    Stage 2: last-resort LLM call to extract and repair a broken JSON array.
    Only triggered when all regex-based repairs have failed.
    Capped at 3000 chars to keep the call cheap and fast.
    """
    print(f"[repair] Calling {repair_llm_name} to repair malformed JSON...")
    llm = get_llm(repair_llm_name)

    prompt = (
        "The following text contains a JSON array of gene objects but it is malformed "
        "or mixed with other text.\n"
        "Extract ONLY the JSON array. Each object must have exactly these keys: "
        '"Gene", "Source Reference", "PMID", "Journal".\n'
        "Return ONLY the valid JSON array, nothing else. "
        "If no gene objects can be found, return an empty array: []\n\n"
        f"Text to repair:\n{raw[:3000]}"
    )

    try:
        result = llm.invoke([
            SystemMessage(content="You extract and repair JSON arrays. Return only valid JSON, no explanation."),
            HumanMessage(content=prompt)
        ])
        return json.loads(result.content.strip())
    except Exception as e:
        print(f"[repair] LLM repair also failed: {e}")
        return []


def safe_json_loads(raw_output: str, repair_llm_name: str = "llama3.1:8b") -> list:
    """
    Parse LLM output into a list of gene objects using a two-stage repair strategy.
    Stage 1: regex (strips <think>, fences, finds array, fixes truncation)
    Stage 2: LLM repair call (last resort)
    """
    # Direct parse — happy path
    try:
        return json.loads(raw_output)
    except json.JSONDecodeError:
        pass

    # Stage 1 — regex repair
    cleaned = repair_json_with_regex(raw_output)
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError:
        pass

    # Stage 1b — fix truncation on the cleaned version
    partial = cleaned
    if partial.count("[") > partial.count("]"):
        partial += "]"
    if partial.count("{") > partial.count("}"):
        partial += "}"
    try:
        return json.loads(partial)
    except json.JSONDecodeError:
        pass

    # Stage 2 — LLM repair (last resort)
    return repair_json_with_llm(raw_output, repair_llm_name)


def build_full_text_context(documents_fulltext: list) -> str:
    """Format full text documents into a context string for the LLM."""
    parts = []
    for d in documents_fulltext:
        sections_text = "\n".join([
            f"  [{s.get('title', section_key)}]: {s.get('text', '')}"
            for section_key, s in d.get("sections", {}).items()
        ])
        parts.append(
            f"PMID: {d.get('pmid')}\n"
            f"Title: {d.get('title')}\n"
            f"Journal: {d.get('journal')}\n"
            f"Full Text:\n{sections_text}"
        )
    return "\n\n".join(parts)


def build_abstract_context(documents: list) -> str:
    """Format abstract documents into a context string for the LLM."""
    return "\n\n".join([
        f"PMID: {d.get('pmid')}\nTitle: {d.get('title')}\n"
        f"Journal: {d.get('journal')}\nAbstract: {d.get('abstract')}"
        for d in documents
    ])


# ── nodes ─────────────────────────────────────────────────────────────────────

def retrieve_pubtator_abstracts(state: LLMGraphState):
    """
    Retrieves abstracts (up to 250) and full texts (up to 50) from PubTator.
    - PMIDs that appear in both are kept only in full texts (more content).
    - Abstracts cache: abstracts/gene_annotated_abstracts/{name}.json
    - Full texts cache: abstracts/full_texts/{name}.json
    Returns: documents (abstracts, overlap removed) + documents_fulltext
    """
    geneset = state["geneset"]
    name = geneset["raw_go_name"].replace("/", " or ")
    print(f"Retrieving documents for gene set: {name}")

    abs_cache  = f"abstracts/gene_annotated_abstracts/{name}.json"
    full_cache = f"abstracts/gene_annotated_full_text/{name}.json"

    # ── load or fetch abstracts ───────────────────────────────────────────────
    if os.path.exists(abs_cache):
        print(f"Cached abstracts found. Loading...")
        with open(abs_cache, "r") as f:
            abstracts = json.load(f)
    else:
        print(f"Downloading abstracts...")
        pmids = Pubtator.search_pubtator_ID(query=name, limit=25)

        abstracts = []
        DELAY = 1.0 / 3
        for pmid in pmids:
            try:
                abs_data = Pubtator.export_abstract(pmid, check_for_genes=False)
                if abs_data:
                    abstracts.append(abs_data)
                time.sleep(DELAY)
            except Exception:
                time.sleep(DELAY)
                continue

        os.makedirs("abstracts/gene_annotated_abstracts", exist_ok=True)
        with open(abs_cache, "w") as f:
            json.dump(abstracts, f, indent=2)
        print(f"Saved {len(abstracts)} abstracts for {name}")

    # ── load or fetch full texts ──────────────────────────────────────────────
    if os.path.exists(full_cache):
        print(f"Cached full texts found. Loading...")
        with open(full_cache, "r") as f:
            full_texts = json.load(f)
    else:
        print(f"Downloading full texts (top {FULL_TEXT_LIMIT})...")
        # Reuse the same PMID order from abstracts, take top N
        full_pmids = [a["pmid"] for a in abstracts[:FULL_TEXT_LIMIT] if a.get("pmid")]

        full_texts = []
        DELAY = 1.0 / 3
        for pmid in full_pmids:
            try:
                ft_data = Pubtator.export_full_text(pmid, check_for_genes=False)
                if ft_data:
                    full_texts.append(ft_data)
                time.sleep(DELAY)
            except Exception:
                time.sleep(DELAY)
                continue

        os.makedirs("abstracts/full_texts", exist_ok=True)
        with open(full_cache, "w") as f:
            json.dump(full_texts, f, indent=2)
        print(f"Saved {len(full_texts)} full texts for {name}")

    # ── deduplicate: remove from abstracts any PMID already in full texts ─────
    full_text_pmids = {str(ft.get("pmid")) for ft in full_texts}
    abstracts_deduped = [
        a for a in abstracts
        if str(a.get("pmid")) not in full_text_pmids
    ]
    removed = len(abstracts) - len(abstracts_deduped)
    if removed:
        print(f"Removed {removed} abstracts that overlap with full texts")

    print(f"Final: {len(abstracts_deduped)} abstracts + {len(full_texts)} full texts")

    return {
        "documents":          abstracts_deduped,
        "documents_fulltext": full_texts
    }


def grade_abstracts(state: LLMGraphState):
    """
    Grades abstracts and saves them to state['documents_filtered'][llm_name].
    Checks filtered_pmids cache first — if PMIDs already graded for this gene
    set, reconstructs filtered docs from the cached list and skips re-grading.
    """
    llm_name = state.get("llm_name", "deepseek-r1:8b")
    print(f"---CHECK ABSTRACT RELEVANCE ({llm_name})---")

    geneset   = state["geneset"]
    gs_name   = geneset["raw_go_name"]
    documents = state.get("documents", [])

    if not documents:
        return {"documents_filtered": {llm_name: []}}

    filtered_pmids_file = f"out/filtered_pmids_{llm_name}.json"

    # ── cache check ──────────────────────────────────────────────────────────
    cache: dict = {}
    if os.path.exists(filtered_pmids_file):
        try:
            with open(filtered_pmids_file, "r") as f:
                cache = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            print(f"[{llm_name}] Could not read abstract PMIDs cache: {e} — re-grading")
            cache = {}

    if gs_name in cache:
        cached_pmids = set(str(p) for p in cache[gs_name])
        filtered = [d for d in documents if str(d.get("pmid")) in cached_pmids]
        print(f"[{llm_name}] Cache hit for '{gs_name}' — "
              f"restored {len(filtered)} filtered abstract(s) from {len(cached_pmids)} cached PMIDs")
        return {"documents_filtered": {llm_name: filtered}}

    # ── grade ────────────────────────────────────────────────────────────────
    question = build_grader_question(geneset)
    llm      = get_llm_json_mode(llm_name)
    filtered = []

    for doc in documents:
        abstract_text = (
            f"Title: {doc.get('title', '')}\n"
            f"Journal: {doc.get('journal', '')}\n"
            f"Abstract: {doc.get('abstract', '')}"
        )
        try:
            result = llm.invoke([
                SystemMessage(content=grade_abstracts_instructions),
                HumanMessage(content=abstract_grader_prompt.format(
                    question=question,
                    document=abstract_text
                ))
            ])
            grade_data = json.loads(result.content)
            if grade_data.get("binary_score", "no").strip().lower() == "yes":
                filtered.append(doc)
        except Exception:
            continue

    print(f"[{llm_name}] Filtered {len(documents)} → {len(filtered)} abstracts")

    # ── persist cache ─────────────────────────────────────────────────────────
    filtered_pmids = [str(d.get("pmid")) for d in filtered if d.get("pmid")]
    os.makedirs(os.path.dirname(filtered_pmids_file) or ".", exist_ok=True)

    existing_cache: dict = {}
    if os.path.exists(filtered_pmids_file):
        try:
            with open(filtered_pmids_file, "r") as f:
                existing_cache = json.load(f)
        except (json.JSONDecodeError, OSError):
            existing_cache = {}

    existing_cache[gs_name] = filtered_pmids
    with open(filtered_pmids_file, "w") as f:
        json.dump(existing_cache, f, indent=2)

    print(f"[{llm_name}] Saved {len(filtered_pmids)} filtered abstract PMIDs "
          f"for '{gs_name}' → {filtered_pmids_file}")
    return {"documents_filtered": {llm_name: filtered}}
    


def grade_full_texts(state: LLMGraphState):
    """
    Grades full texts for gene set relevance, one document at a time.
    Checks filtered_pmids_fulltext cache first — if PMIDs already graded for
    this gene set, restores filtered docs from the cache and skips re-grading.
    Saves filtered full texts to state['documents_fulltext_filtered'][llm_name].
    """
    llm_name = state.get("llm_name", "deepseek-r1:8b")
    print(f"---CHECK FULL TEXT RELEVANCE ({llm_name})---")

    geneset            = state["geneset"]
    gs_name            = geneset["raw_go_name"]
    documents_fulltext = state.get("documents_fulltext", [])

    if not documents_fulltext:
        print(f"[{llm_name}] No full texts to grade for '{gs_name}'")
        return {"documents_fulltext_filtered": {llm_name: []}}

    filtered_pmids_file = f"out/filtered_pmids_fulltext_{llm_name}.json"

    # ── cache check ──────────────────────────────────────────────────────────
    cache: dict = {}
    if os.path.exists(filtered_pmids_file):
        try:
            with open(filtered_pmids_file, "r") as f:
                cache = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            print(f"[{llm_name}] Could not read full text PMIDs cache: {e} — re-grading")
            cache = {}

    if gs_name in cache:
        cached_pmids = set(str(p) for p in cache[gs_name])
        filtered = [d for d in documents_fulltext if str(d.get("pmid")) in cached_pmids]
        print(f"[{llm_name}] Cache hit for '{gs_name}' — "
              f"restored {len(filtered)} filtered full text(s) from {len(cached_pmids)} cached PMIDs")
        return {"documents_fulltext_filtered": {llm_name: filtered}}

    # ── grade each full text individually ────────────────────────────────────
    question = build_grader_question(geneset)
    llm      = get_llm_json_mode(llm_name)
    filtered = []

    for doc in documents_fulltext:
        doc_pmid = doc.get("pmid", "?")

        # Flatten sections into a single context string
        sections_text = "\n\n".join([
            f"Section: {s.get('title', k)}\n{s.get('text', '')}"
            for k, s in doc.get("sections", {}).items()
        ])
        full_text_body = (
            f"Title: {doc.get('title', '')}\n"
            f"Journal: {doc.get('journal', '')}\n\n"
            f"{sections_text}"
        )

        try:
            result = llm.invoke([
                SystemMessage(content=grade_full_texts_instructions),
                HumanMessage(content=abstract_grader_prompt.format(
                    question=question,
                    document=full_text_body
                ))
            ])
            grade_data = json.loads(result.content)
            if grade_data.get("binary_score", "no").strip().lower() == "yes":
                filtered.append(doc)
                print(f"[{llm_name}] PMID {doc_pmid} ✓ relevant")
            else:
                print(f"[{llm_name}] PMID {doc_pmid} ✗ not relevant")
        except Exception as e:
            print(f"[{llm_name}] PMID {doc_pmid} grading error: {e} — skipping")
            continue

    print(f"[{llm_name}] Filtered {len(documents_fulltext)} → {len(filtered)} full texts")

    # ── persist cache ─────────────────────────────────────────────────────────
    filtered_pmids = [str(d.get("pmid")) for d in filtered if d.get("pmid")]
    os.makedirs(os.path.dirname(filtered_pmids_file) or ".", exist_ok=True)

    existing_cache: dict = {}
    if os.path.exists(filtered_pmids_file):
        try:
            with open(filtered_pmids_file, "r") as f:
                existing_cache = json.load(f)
        except (json.JSONDecodeError, OSError):
            existing_cache = {}

    existing_cache[gs_name] = filtered_pmids
    with open(filtered_pmids_file, "w") as f:
        json.dump(existing_cache, f, indent=2)

    print(f"[{llm_name}] Saved {len(filtered_pmids)} filtered full text PMIDs "
          f"for '{gs_name}' → {filtered_pmids_file}")
    return {"documents_fulltext_filtered": {llm_name: filtered}}


FULL_TEXT_CHUNK_SIZE = 2   # full texts per LLM call — they're large
ABSTRACT_CHUNK_SIZE  = 10   # abstracts per LLM call


def _run_chunk(
    chunk_context: str,
    question: str,
    llm_name: str,
    llm,
    json_outfile: str,
    raw_outfile: str,
    repair_llm: str = "llama3.1:8b",
) -> list:
    """
    Send one context chunk to the LLM, parse the result, and append-safe-save.
    Returns the parsed gene list for this chunk (may be empty).
    Only writes to the raw file on a genuine parse failure — not on [] responses.
    """
    messages = [
        SystemMessage(content=(
            "You are a precise biomedical text mining assistant specialized in "
            "Triticum aestivum (common wheat) gene set curation. "
            "Respond only in valid JSON. "
            "Only report genes that are EXPLICITLY named in the provided text."
        )),
        HumanMessage(content=rag_prompt.format(context=chunk_context, question=question))
    ]

    try:
        result     = llm.invoke(messages)
        raw_output = result.content.strip()

        # Strip <think>...</think> blocks and markdown fences before parsing
        cleaned_output = re.sub(r"<think>.*?</think>", "", raw_output, flags=re.DOTALL).strip()
        cleaned_output = re.sub(r"```(?:json)?", "", cleaned_output).strip(" `\n")

        # Explicit empty array — model found no genes, valid response, don't log as error
        if cleaned_output == "[]":
            print(f"[{llm_name}] Chunk: no wheat genes found (explicit [])")
            return []

        generation = safe_json_loads(cleaned_output, repair_llm_name=repair_llm)

        if not generation:
            generation = safe_json_loads(
                clean_model_output(cleaned_output), repair_llm_name=repair_llm
            )

        if not generation:
            # Genuine parse failure — save raw output for inspection
            print(f"[{llm_name}] Chunk: could not parse output — saved to raw file")
            with open(raw_outfile, "a") as f:
                f.write(f"\n\n--- PARSE FAILURE ---\n{raw_output}")
            return []

    except Exception as e:
        print(f"[{llm_name}] Chunk invocation error: {e}")
        with open(raw_outfile, "a") as f:
            f.write(f"\n\n--- ERROR ---\n{e}")
        return []

    added = append_safe_save(generation, json_outfile)
    print(f"[{llm_name}] Chunk: {len(generation)} gene(s) found, {added} new")
    return generation


def generate(state: LLMGraphState):
    """
    Generates gene list by processing documents in small chunks to stay within
    context limits and prevent hallucination from information overload.

    Full texts: FULL_TEXT_CHUNK_SIZE per call (large content)
    Abstracts:  ABSTRACT_CHUNK_SIZE per call

    Uses append_safe_save to deduplicate across all chunks on (Gene, PMID).
    """
    llm_name = state.get("llm_name", "deepseek-r1:8b")
    print(f"---GENERATE ({llm_name})---")

    geneset   = state["geneset"]
    safe_name = geneset["raw_go_name"].replace("/", " or ")

    all_filtered       = state.get("documents_filtered", {})
    documents          = all_filtered.get(llm_name, [])
    all_ft_filtered    = state.get("documents_fulltext_filtered", {})
    documents_fulltext = all_ft_filtered.get(llm_name, [])

    if not documents and not documents_fulltext:
        print(f"[{llm_name}] No documents available")
        return {"generation": {llm_name: []}}

    question = (
        f"Identify genes explicitly mentioned in the provided text that are "
        f"associated with the Triticum aestivum (common wheat) gene set '{safe_name}'. "
        f"Definition: '{geneset['go_term_definition']}'."
    )

    out_dir = f"out/geneset_generations/{llm_name}"
    os.makedirs(out_dir, exist_ok=True)
    json_outfile = f"{out_dir}/{safe_name}.json"
    raw_outfile  = f"{out_dir}/{safe_name}_raw.txt"

    # Initialise output files fresh for this run
    with open(raw_outfile, "w") as f:
        f.write(f"Raw LLM outputs for gene set: {geneset['raw_go_name']} ({llm_name})\n")
    with open(json_outfile, "w") as f:
        json.dump([], f, indent=2)

    llm        = get_llm(llm_name)
    repair_llm = "llama3.1:8b"
    all_genes  = []

    # ── full text chunks ──────────────────────────────────────────────────────
    for i in range(0, len(documents_fulltext), FULL_TEXT_CHUNK_SIZE):
        chunk = documents_fulltext[i : i + FULL_TEXT_CHUNK_SIZE]
        print(f"[{llm_name}] Full text chunk {i // FULL_TEXT_CHUNK_SIZE + 1} "
              f"({len(chunk)} doc(s), PMIDs: {[d.get('pmid') for d in chunk]})")
        context = "=== FULL TEXT SOURCES ===\n" + build_full_text_context(chunk)
        genes   = _run_chunk(context, question, llm_name, llm, json_outfile, raw_outfile, repair_llm)
        all_genes.extend(genes)

    # ── abstract chunks ───────────────────────────────────────────────────────
    for i in range(0, len(documents), ABSTRACT_CHUNK_SIZE):
        chunk = documents[i : i + ABSTRACT_CHUNK_SIZE]
        print(f"[{llm_name}] Abstract chunk {i // ABSTRACT_CHUNK_SIZE + 1} "
              f"({len(chunk)} doc(s), PMIDs: {[d.get('pmid') for d in chunk]})")
        context = "=== ABSTRACT SOURCES ===\n" + build_abstract_context(chunk)
        genes   = _run_chunk(context, question, llm_name, llm, json_outfile, raw_outfile, repair_llm)
        all_genes.extend(genes)

    print(f"[{llm_name}] Generation complete — {len(all_genes)} total gene(s) across all chunks")
    return {"generation": {llm_name: all_genes}}

def verify_metadata(state: LLMGraphState):
    """
    Qwen post-processing node: corrects PMID and Journal fields in the generation
    by cross-referencing against the actual documents (full texts + abstracts).
    If a match cannot be found, the original value is left unchanged.
    Only modifies PMID and Journal — gene names and source references are untouched.
    """
    print(f"---VERIFY METADATA (qwen)---")

    llm_name  = state.get("llm_name", "deepseek-r1:8b")
    geneset   = state["geneset"]
    safe_name = geneset["raw_go_name"].replace("/", " or ")

    # Get the current generation for this LLM
    all_generation     = state.get("generation", {})
    generation         = all_generation.get(llm_name, [])
    documents          = state.get("documents", [])
    documents_fulltext = state.get("documents_fulltext", [])

    if not generation:
        print("[verify] No generation to verify.")
        return {"generation": {llm_name: []}}

    # Build a compact document reference list for Qwen (PMID + title + journal only)
    # Keeps the prompt short — Qwen only needs to match, not read full content
    ref_lines = []
    for d in documents_fulltext + documents:
        pmid    = d.get("pmid", "")
        title   = d.get("title", "")
        journal = d.get("journal", "")
        if pmid:
            ref_lines.append(f'PMID: {pmid} | Journal: {journal} | Title: {title}')
    references = "\n".join(ref_lines)

    prompt = (
        "You are a metadata verification assistant. Below is a list of gene objects "
        "extracted from scientific literature, followed by the actual source documents.\n\n"
        "Your task: for each gene object, check whether the PMID and Journal fields "
        "match a real document in the source list.\n"
        "- If the PMID matches a document, ensure Journal is also correct.\n"
        "- If the PMID cannot be found in the source list, leave the original PMID and Journal unchanged.\n"
        "- Do NOT change Gene or Source Reference fields.\n"
        "- Return the full corrected list as a valid JSON array with the same structure.\n\n"
        f"SOURCE DOCUMENTS:\n{references}\n\n"
        f"GENE OBJECTS TO VERIFY:\n{json.dumps(generation, indent=2)}"
    )

    llm = get_llm(VERIFY_MODEL)
    verified = generation  # fallback to original if Qwen fails

    try:
        result = llm.invoke([
            SystemMessage(content=(
                "You verify and correct metadata in JSON arrays. "
                "Return only a valid JSON array, no explanation."
            )),
            HumanMessage(content=prompt)
        ])
        raw       = result.content.strip()
        parsed    = safe_json_loads(raw, repair_llm_name="llama3.1:8b")
        if parsed:
            verified = parsed
            print(f"[verify] Qwen verified {len(verified)} gene object(s)")
        else:
            print("[verify] Qwen output could not be parsed — keeping original metadata")
    except Exception as e:
        print(f"[verify] Qwen error: {e} — keeping original metadata")
    
    print(f"[verify] Final verified generation:\n{json.dumps(verified, indent=2)}")

    # Overwrite the JSON file with verified metadata
    json_outfile = f"out/geneset_generations/{llm_name}/{safe_name}.json"
    if os.path.exists(json_outfile) and verified:
        # Reload existing, replace entries that match by Gene+PMID key
        try:
            with open(json_outfile, "r") as f:
                existing = json.load(f)

            # Build lookup of verified objects by (Gene, original PMID)
            verified_lookup = {
                str(g.get("Gene", "")).strip().lower(): g
                for g in verified
            }
            # Update existing entries with corrected metadata
            for entry in existing:
                print(f"\n\nGot here\n\n")
                gene_key = str(entry.get("Gene", "")).strip().lower()
                if gene_key in verified_lookup:
                    entry["PMID"]    = verified_lookup[gene_key].get("PMID", entry["PMID"])
                    entry["Journal"] = verified_lookup[gene_key].get("Journal", entry["Journal"])

            with open(json_outfile, "w") as f:
                json.dump(existing, f, indent=2)
            print(f"[verify] Updated metadata saved to {json_outfile}")
        except Exception as e:
            print(f"[verify] Could not update file: {e}")

    return {"generation": {llm_name: verified}}


def create_control_flow():
    workflow = StateGraph(LLMGraphState)

    workflow.add_node("retrieve",          retrieve_pubtator_abstracts)
    workflow.add_node("grade",             lambda s: grade_abstracts(s))
    workflow.add_node("grade_full_texts",  lambda s: grade_full_texts(s))
    workflow.add_node("generate",          lambda s: generate(s))
    workflow.add_node("verify_metadata",   lambda s: verify_metadata(s))

    workflow.set_entry_point("retrieve")
    workflow.add_edge("retrieve",         "grade")
    workflow.add_edge("grade",            "grade_full_texts")
    workflow.add_edge("grade_full_texts", "generate")
    workflow.add_edge("generate",          END)

    return workflow.compile()