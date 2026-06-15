"""
convert_genes.py
----------------
Pre-BLAST gene conversion pipeline.

For each gene extracted by the maker LLM:
  1. Try BioMart chip table lookup (instant)
  2. If not found: small LLM normalises name, retry chip table
  3. If still not found: LLM classifies is_gene (yes/no)
     - is_gene=True  → store, tag as 'unconverted_gene' (queued for BLAST)
     - is_gene=False → store, tag as 'not_gene'
  4. Write conversion_results.json + lightweight HTML viewer

Usage:
    python convert_genes.py --llm deepseek-r1:8b
    python convert_genes.py --llm qwen3:32b --biomart data/merged_biomart.tsv
"""

import os
import re
import json
import csv
import argparse
from pathlib import Path
from collections import defaultdict
from langchain_core.messages import HumanMessage, SystemMessage

# ── config ────────────────────────────────────────────────────────────────────
DEFAULT_BIOMART   = "data/merged_biomart.tsv"
DEFAULT_GEN_BASE  = "out/geneset_generations"
OUTPUT_DIR        = "out/conversion"
RESULTS_JSON      = "out/conversion/conversion_results.json"
VIEWER_HTML       = "out/conversion/viewer.html"
NORMALIZE_LLM     = "deepseek-r1:8b"
MAKER_FILE       = "out/processed_genesets_maker_deepseek-r1:8b.txt"

# Tags
TAG_CONVERTED           = "converted"
TAG_CONVERTED_NORMALISED = "converted_normalised"
TAG_UNCONVERTED_GENE    = "unconverted_gene"
TAG_NOT_GENE            = "not_gene"
# ─────────────────────────────────────────────────────────────────────────────


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
        "transcript_id":     "Transcript stable ID",
        "gene_name":         "Gene name",
        "gene_synonym":      "Gene Synonym",
        "ncbi_id":           "NCBI gene (formerly Entrezgene) ID",
        "ncbi_accession":    "NCBI gene (formerly Entrezgene) accession",
        "refseq_mrna":       "RefSeq mRNA ID",
        "uniprot_swiss":     "UniProtKB/Swiss-Prot ID",
        "uniprot_trembl":    "UniProtKB/TrEMBL ID",
        "uniprot_symbol":    "UniProtKB Gene Name symbol",
        "wikigene":          "WikiGene name",
        "interpro_desc":     "InterPro Short Description",
        "protein_domain":    "Protein domain description",
        "ncbi_gene_desc":    "NCBI gene (formerly Entrezgene) description",
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


def chip_lookup(name: str, lookups: dict) -> tuple[str | None, str]:
    """
    Try all lookup tables. Returns (TraesCS_ID, source) or (None, 'not_found').
    """
    name_lower = name.strip().lower()
    name_exact = name.strip()

    # already in TraesCS format
    if re.match(r"^TraesCS\w+", name_exact):
        return name_exact.split(".")[0], "direct_traescs"

    # old Traes_ MIPS format — NOT a proper TraesCS ID.
    # Do NOT short-circuit; let it fall through to BioMart tables.
    # If BioMart has no mapping it will come back not_found → LLM → unconverted_gene.

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
        table = lookups.get(table_key, {})
        result = table.get(name_exact) or table.get(name_lower)
        if result:
            return result, source_label

    return None, "not_found"


# ── LLM normalise + classify ──────────────────────────────────────────────────

def llm_normalise_and_classify(
    gene_name: str,
    llm_name:  str = NORMALIZE_LLM
) -> dict:
    """
    Single LLM call that:
    1. Normalises the gene name (strips species prefixes, fixes case etc.)
    2. Suggests alternative search terms
    3. Determines whether this is actually a gene

    Returns dict with: normalised_name, alternatives, is_gene, confidence
    """
    try:
        from utils import get_llm_json_mode
        llm = get_llm_json_mode(llm_name)
    except Exception as e:
        print(f"  [LLM] Could not load model {llm_name}: {e}")
        return {
            "normalised_name": gene_name,
            "alternatives":    [],
            "is_gene":         True,
            "confidence":      "low"
        }

    prompt = f"""Analyse this identifier from a wheat genomics abstract.

Identifier: "{gene_name}"

Return ONLY this JSON:
{{
  "normalised_name": "cleaned version for database lookup (strip Ta prefix, standardise case, remove subgenome suffix)",
  "alternatives": ["list", "of", "other", "names", "to", "try"],
  "is_gene": true or false,
  "confidence": "high" or "medium" or "low"
}}

Rules for is_gene:
- true if it is a gene, protein, enzyme, or RNA
- false if it is a GO term, phenotype description, tissue name, chemical compound, or biological process

Rules for alternatives (include up to 4):
- strip Ta/Tt species prefix e.g. TaHSP70 → HSP70
- strip subgenome suffix e.g. HSP70-B1 → HSP70
- expand abbreviations if known e.g. VRN1 → VERNALIZATION1
- include Arabidopsis ortholog name if known"""

    try:
        result = llm.invoke([
            SystemMessage(content="You are a wheat genomics expert. Return only valid JSON."),
            HumanMessage(content=prompt)
        ])
        raw = re.sub(r"<think>.*?</think>", "", result.content, flags=re.DOTALL).strip()
        # extract JSON from fences if present
        match = re.search(r"```(?:json)?\s*([\s\S]*?)```", raw)
        if match:
            raw = match.group(1).strip()
        return json.loads(raw)
    except Exception as e:
        print(f"  [LLM] Parse error: {e}")
        return {
            "normalised_name": gene_name,
            "alternatives":    [],
            "is_gene":         True,
            "confidence":      "low"
        }


# ── per-gene conversion ───────────────────────────────────────────────────────

def convert_gene(
    gene_obj:  dict,
    lookups:   dict,
    llm_name:  str
) -> dict:
    """
    Convert a single gene object. Returns enriched dict with conversion metadata.
    """
    gene_name = str(gene_obj.get("Gene", "")).strip()

    result = {
        "gene_input":        gene_name,
        "traescs_id":        None,
        "conversion_source": None,
        "tag":               None,
        "is_gene":           None,
        "normalised_name":   None,
        "llm_used":          False,
        "pmid":              str(gene_obj.get("PMID", "")).strip(),
        "journal":           str(gene_obj.get("Journal", "")).strip(),
        "source_reference":  str(gene_obj.get("Source Reference", "")).strip(),
    }

    if not gene_name:
        result["tag"] = TAG_NOT_GENE
        return result

    # Stage 1 — direct chip lookup
    traescs, source = chip_lookup(gene_name, lookups)
    if traescs:
        result.update({
            "traescs_id":        traescs,
            "conversion_source": source,
            "tag":               TAG_CONVERTED,
            "is_gene":           True
        })
        return result

    # Stage 2 — LLM normalise + retry chip
    print(f"  [normalise] {gene_name}")
    classification  = llm_normalise_and_classify(gene_name, llm_name)
    normalised      = classification.get("normalised_name", gene_name)
    alternatives    = classification.get("alternatives", [])
    is_gene         = classification.get("is_gene", True)

    result["normalised_name"] = normalised
    result["llm_used"]        = True
    result["is_gene"]         = is_gene

    # retry chip with normalised name and alternatives
    for candidate in [normalised] + alternatives:
        traescs, source = chip_lookup(candidate, lookups)
        if traescs:
            result.update({
                "traescs_id":        traescs,
                "conversion_source": source,
                "tag":               TAG_CONVERTED_NORMALISED,
            })
            return result

    # Stage 3 — classify and tag
    if not is_gene:
        result["tag"] = TAG_NOT_GENE
    else:
        result["tag"] = TAG_UNCONVERTED_GENE   # queued for BLAST later

    return result


# ── main pipeline ─────────────────────────────────────────────────────────────

def run(llm_name: str, biomart_file: str):
    gen_dir = Path(DEFAULT_GEN_BASE) / llm_name
    if not gen_dir.exists():
        print(f"Generation directory not found: {gen_dir}")
        return

    json_files = sorted(gen_dir.glob("*.json"))
    # only use json files in the maker file (i.e. those that have been processed/verified)
    if Path(MAKER_FILE).exists():
        with open(MAKER_FILE, "r") as f:
            processed_genesets = {line.strip() for line in f if line.strip()}
        json_files = [f for f in json_files if f.stem in processed_genesets]
    print(f"\nFound {len(json_files)} gene set files in {gen_dir}\n")

    print("Building lookup tables...")
    lookups = build_lookup_tables(biomart_file)
    print()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_results  = []
    stats = {
        TAG_CONVERTED:            0,
        TAG_CONVERTED_NORMALISED: 0,
        TAG_UNCONVERTED_GENE:     0,
        TAG_NOT_GENE:             0,
        "total_genes":            0,
        "total_gene_sets":        0
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

        for gene_obj in genes:
            result = convert_gene(gene_obj, lookups, llm_name)
            result["gene_set"] = gene_set_name
            geneset_results.append(result)
            stats["total_genes"] += 1
            tag = result.get("tag", TAG_NOT_GENE)
            if tag in stats:
                stats[tag] += 1

        all_results.extend(geneset_results)

        # save incrementally every 50 gene sets
        # save incrementally every 50 gene sets
        if i % 50 == 0:
            partial = [
                r for r in all_results
                if r.get("gene_input", "").strip() and r.get("tag") != TAG_NOT_GENE
            ]
            with open(RESULTS_JSON, "w") as f:
                json.dump(partial, f, indent=2)

    # final save
    # final save — strip not_gene and empty gene_input entries
    filtered_results = [
        r for r in all_results
        if r.get("gene_input", "").strip()          # drop blank gene_input
        and r.get("tag") != TAG_NOT_GENE            # drop LLM-confirmed non-genes
    ]
    with open(RESULTS_JSON, "w") as f:
        json.dump(filtered_results, f, indent=2)

    # generate HTML viewer
    generate_viewer(RESULTS_JSON, VIEWER_HTML)

    print(f"\n{'─' * 50}")
    print(f"SUMMARY")
    print(f"  Gene sets processed   : {stats['total_gene_sets']}")
    print(f"  Total genes           : {stats['total_genes']}")
    print(f"  Converted (direct)    : {stats[TAG_CONVERTED]}")
    print(f"  Converted (normalised): {stats[TAG_CONVERTED_NORMALISED]}")
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
  .tag-converted_normalised { background: #1c3a14; color: #86efac; }
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
    <option value="converted_normalised">Converted ✓ (normalised)</option>
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
        <th>Normalised</th>
        <th>TraesCS ID</th>
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

// Only two tag flags: converted and gene ⚑
// not_gene and empty entries are stripped from the JSON before this viewer is generated
const TAG_LABELS = {
  "converted":            { label: "converted ✓", cls: "tag-converted"  },
  "converted_normalised": { label: "converted ✓", cls: "tag-converted"  },
  "unconverted_gene":     { label: "gene ⚑",       cls: "tag-gene"       },
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
    converted:            ['Converted',       'green'],
    converted_normalised: ['Normalised',      'green'],
    unconverted_gene:     ['Queued for BLAST','amber'],
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
    const matchQ   = !q || [r.gene_input, r.normalised_name, r.traescs_id,
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
      <td>${r.normalised_name && r.normalised_name !== r.gene_input ? esc(r.normalised_name) : '<span style="color:var(--muted)">—</span>'}</td>
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
        description="Pre-BLAST gene conversion: BioMart lookup + LLM normalisation"
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
        "--normalize-llm",
        default=NORMALIZE_LLM,
        help=f"LLM to use for name normalisation (default: {NORMALIZE_LLM})"
    )
    args = parser.parse_args()

    run(
        llm_name     = args.llm,
        biomart_file = args.biomart
    )