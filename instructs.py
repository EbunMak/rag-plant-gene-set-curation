### LLM
from langchain_ollama import ChatOllama

local_llm = "llama3.1:8b"
llm = ChatOllama(model=local_llm, temperature=0)
llm_json_mode = ChatOllama(model=local_llm, temperature=0, format="json")


# -------------------------------------------------------------------
# GRADER 1 — Gene-set/phenotype relevance filter (used in rag_pipeline)
# -------------------------------------------------------------------
# Purpose: keep abstracts that mention the GO term / gene set topic AND
#          have any connection to wheat (T. aestivum) or closely related
#          Triticeae / Poaceae species.
#
# Key changes from original:
#   - Separates the two acceptance criteria clearly (topic + wheat link)
#   - Provides explicit PASS examples so a small LLM can calibrate
#   - Removes the contradictory hedge "even if not explicitly mentioned"
#   - Instructs the model to be INCLUSIVE rather than conservative,
#     because false-negatives are more costly than false-positives here
#     (a second grader or the generate step can discard irrelevant genes)
# -------------------------------------------------------------------

grade_abstracts_instructions = """You are a relevance filter for a wheat gene-set curation pipeline.

Your job is to decide whether a scientific abstract should be KEPT or DISCARDED.

KEEP the abstract (binary_score: "yes") if it satisfies ALL of the following:
  1. TOPIC MATCH — The abstract discusses a biological process, molecular function,
     cellular component, phenotype, or pathway that is related to the gene set topic
     provided in the question (e.g. "cell wall organisation", "response to drought").
  2. WHEAT LINK — The abstract involves Triticum aestivum, T. durum, T. dicoccoides,
     Aegilops tauschii, or another Triticeae / Poaceae species (barley, rice, maize,
     Arabidopsis are acceptable proxies when functional conservation with wheat
     is plausible).

DISCARD the abstract (binary_score: "no") ONLY if:
  - The topic is clearly unrelated to the gene set (e.g. a cancer drug trial when
    the gene set is "starch biosynthesis"), OR
  - The abstract contains no plant biology whatsoever (e.g. purely human/animal study
    with no transferable functional context).

IMPORTANT — Be INCLUSIVE. When in doubt, return "yes".
A downstream step will verify individual gene-phenotype claims.
Rejecting a borderline abstract loses evidence permanently.

Respond ONLY with valid JSON: {"binary_score": "yes"} or {"binary_score": "no"}.
Do not add any explanation or extra keys."""


# -------------------------------------------------------------------
# GRADER 2 — Gene + phenotype co-mention filter (validation pipeline)
# -------------------------------------------------------------------
# Unchanged structurally, but tightened language for small LLMs.
# -------------------------------------------------------------------

grade_abstracts_instructions2 = """You are a relevance filter for a wheat gene validation pipeline.

Decide whether the abstract discusses BOTH:
  1. The GENE named in the question (by name, symbol, or clear synonym).
  2. The PHENOTYPE or biological process named in the question.

Return "yes" only when both are discussed in a functionally related context
(same experiment, pathway, mechanism, or association).
Return "no" if either the gene or the phenotype is absent or only mentioned
incidentally with no functional link.

Respond ONLY with valid JSON: {"binary_score": "yes"} or {"binary_score": "no"}.
Do not add any explanation or extra keys."""


# GRADER for full text sections in the gene checker pipeline — same criteria as above but applied to full text instead of abstracts
grade_full_texts_instructions = """You are a relevance filter for a wheat gene validation pipeline that uses full text sections from scientific papers.
Decide whether the full text section discusses BOTH:
  1. The GENE named in the question (by name, symbol, or clear synonym).
  2. The PHENOTYPE or biological process named in the question.

Return "yes" only when both are discussed in a functionally related context
(same experiment, pathway, mechanism, or association).
Return "no" if either the gene or the phenotype is absent or only mentioned
incidentally with no functional link.
Respond ONLY with valid JSON: {"binary_score": "yes"} or {"binary_score": "no"}.
Do not add any explanation or extra keys."""

# -------------------------------------------------------------------
# USER-TURN GRADER PROMPT (used when invoking as a HumanMessage)
# -------------------------------------------------------------------
# This is used by the abstract_grader_prompt template variable in the
# pipeline. Kept short — the system message carries the criteria.
# -------------------------------------------------------------------

abstract_grader_prompt = """Gene set topic: {question}

Abstract to evaluate:
{document}

Does this abstract satisfy the KEEP criteria described in your instructions?
Return JSON: {{"binary_score": "yes"}} or {{"binary_score": "no"}}."""


# -------------------------------------------------------------------
# QUESTION STRING — built inside grade_abstracts() in rag_pipeline
# -------------------------------------------------------------------
# Replace the inline question string in grade_abstracts() with this
# helper so the wording is consistent and easy to update in one place.
# -------------------------------------------------------------------

def build_grader_question(geneset: dict) -> str:
    """
    Returns a concise, unambiguous question string for the abstract grader.
    Replaces the convoluted inline string in grade_abstracts().

    Usage in rag_pipeline_gene_set_maker_llm.py:
        from instructs import build_grader_question
        question = build_grader_question(geneset)
    """
    return (
        f"Gene set topic: '{geneset['raw_go_name']}' defined as '{geneset['go_term_definition']}', in Triticum aestivum (wheat) or related Triticeae/Poaceae species. "
        f"Does this abstract contain information about genes, pathways, or molecular functions "
        f"relevant to this gene set?"
    )


# -------------------------------------------------------------------
# RAG GENERATION PROMPT
# -------------------------------------------------------------------

rag_prompt = """You are an assistant for a gene set curation task in Triticum aestivum (common wheat).
 
Task Overview:
Your role is to extract wheat genes that are EXPLICITLY NAMED in the provided text and are relevant to the following biological process:
 
{question}
 
Source text:
{context}
 
CRITICAL RULES — read carefully before responding:
1. ONLY report genes that are explicitly named, by any identifier, in the text above.
2. Do NOT infer, suggest, predict, or hallucinate genes that are not directly stated.
3. If the text contains no Triticum aestivum (wheat) genes, return an empty array: []
4. Every "Source Reference" must be a DIRECT QUOTE from the text — copy the exact sentence(s) that name the gene.
5. Every PMID must come from the document headers in the text above — do not invent or guess PMIDs.
6. ONLY report genes from Triticum aestivum (common wheat). Do NOT report genes from other species
   (e.g. mouse, human, Arabidopsis, rice) even if they are mentioned in the same text. If a gene
   name appears without a species context, only include it if the surrounding text clearly associates
   it with wheat.
Wheat genes may appear under any of the following naming schemes — recognise ALL of them:
 
1. IWGSC gene model IDs (current standard):
   - Gene-level:       TraesCS2B02G154600
   - Transcript-level: TraesCS2B02G154600.1
   Format: TraesCS + chromosome number + subgenome letter (A/B/D) + 02 + G + 6-digit number
 
2. MIPSv2.2 / Phytozome IDs (older assembly):
   - Traes_2BS_1CFB331E9, Traes_2BS_1CFB331E9.1
   Format: Traes_ + chromosome + subgenome + arm (S/L) + underscore + hash
 
3. WGC (Wheat Gene Catalogue) names — phenotype-based:
   - Vrn-A1, Vrn-B1, Vrn-D1 (vernalization)
   - Rht-B1, rht-b1 (reduced height)
   - Lr34, Lr46 (leaf rust resistance)
   - Ppd-D1 (photoperiod)
   Format: trait abbreviation + subgenome letter + number
 
4. Arabidopsis homology-based names:
   - TaAP1-A1, TaAP1-B1, TaAP1-D1, TaSEP3-D1, TaVRN1
   Format: Ta + gene name from model organism + subgenome
 
5. Protein/enzyme family names:
   - TaHSP70-B1, TaMAPK-A3, TaWRKY-D6, TaDREB1
   Format: Ta + protein family abbreviation + optional subgenome + copy number
 
6. Affymetrix probe IDs:
   - TaAffx.98024.1.A1_at
   Format: TaAffx. + numeric ID + version + subgenome + _at
 
7. GenBank accession numbers appearing alongside gene names:
   - AY188331, MK577897
 
If a gene appears under multiple naming schemes in the same document, report it ONCE using the most specific identifier available, in this order of preference:
TraesCS ID > WGC name > homology-based name > protein family name > Affymetrix ID > accession number
 
Output format:
Respond ONLY with a JSON array. Each object must have EXACTLY these four keys:
 
{{
  "Gene": "The gene identifier exactly as it appears in the text, using the most specific naming scheme.",
  "Source Reference": "Direct verbatim quote from the text that names this gene.",
  "PMID": "The PubMed ID from the document header where this gene was found.",
  "Journal": "The journal name from the document header where this gene was found."
}}
 
No duplicate genes. No commentary. No markdown. No extra text. Valid JSON only.
If no wheat genes are explicitly named in the text, return: []
"""
 

rag_prompt2 = """You are an assistant for a gene set validation task in Triticum aestivum (common wheat).  
Here is the context to use to answer the question:
{context}

Task Overview:  
Your role is to determine whether a wheat gene provided in the question is supported by the provided context in relation to a given wheat biological process or disease also in the question.  
Here is the question containing the gene and the biological process or disease:
{question}

Then return JSON with exactly these keys:
1. Gene: gene name  
2. Validation: "yes" or "no" score to indicate whether the gene is validated with the context  
3. Supporting Extract: "Direct quote(s) from the context that validates the gene." and short reasoning (≤5 lines) if more than one abstract is used  
4. PMIDS: list of pmids of the abstract(s) where the inference was made. If the answer for validation is yes then at least one PMID must be provided.  

Do not return any explanation or extra text—only JSON.
"""
