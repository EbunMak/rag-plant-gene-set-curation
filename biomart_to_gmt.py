"""
BioMart TSV → GMT (x3) + JSON converter for wheat gene sets
------------------------------------------------------------
Input:  BioMart TSV with columns:
        Gene stable ID, Transcript stable ID, Gene name, Transcript name,
        Gene Synonym, GO term name, GO term definition, GO term accession,
        GO term evidence code, GO domain

Outputs:
  wheat_gmt_gene_stable_id.gmt       — genes as TraesCS IDs
  wheat_gmt_transcript_stable_id.gmt — genes as transcript IDs (TraesCS...1)
  wheat_gmt_gene_name.gmt            — genes as functional names (WM12, ALMT1)
  go_terms.json                      — metadata for each GO term
"""

import pandas as pd
import json
import sys
from pathlib import Path


# ── config ────────────────────────────────────────────────────────────────────
INPUT_FILE      = "data/mart_exports/mart_export_2.1.tsv"
GMT_GENE_STABLE = "out/wheat_gmt_gene_stable_id_2.1.gmt"
GMT_TRANSCRIPT  = "out/wheat_gmt_transcript_stable_id_2.1.gmt"
GMT_GENE_NAME   = "out/wheat_gmt_gene_name_2.1.gmt"
JSON_OUT        = "out/go_terms_2.1.json"
# ─────────────────────────────────────────────────────────────────────────────


def normalise_go_name(name: str) -> str:
    """'protein targeting to Golgi' → 'PROTEIN_TARGETING_TO_GOLGI'"""
    return name.strip().upper().replace(" ", "_").replace("-", "_")


def load_and_validate(path: str) -> pd.DataFrame:
    sep = "\t" if path.endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    df.columns = df.columns.str.strip()

    required = {
        "Gene stable ID", "Transcript stable ID", "Gene name",
        "GO term name", "GO term definition",
        "GO term accession", "GO term evidence code", "GO domain"
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing columns: {missing}\nFound: {list(df.columns)}"
        )

    before = len(df)
    df = df.dropna(subset=["GO term name", "Gene stable ID"])
    dropped = before - len(df)
    if dropped:
        print(f"   Dropped {dropped} rows with missing GO term or Gene ID.")

    for col in ["Gene stable ID", "Transcript stable ID",
                "Gene name", "GO term name", "GO domain",
                "GO term accession"]:
        df[col] = df[col].str.strip()

    return df


def build_gmt(df: pd.DataFrame, id_col: str) -> dict:
    """
    Build a gene-set dict keyed by normalised GO term name.
    id_col: which column to use as the gene identifier.

    Returns dict:
      normalised_go_name -> {
          "go_domain":  biological_process | molecular_function | cellular_component,
          "go_id":      GO:XXXXXXX,
          "genes":      ordered list of unique IDs from id_col,
          "definition": full GO term definition,
          "evidence":   sorted list of evidence codes,
          "raw_name":   original GO term name
      }
    """
    gmt = {}

    for _, row in df.iterrows():
        go_name  = row["GO term name"]
        gene_id  = row[id_col]
        go_def   = row.get("GO term definition", "")
        go_id    = row.get("GO term accession", "")
        go_dom   = row.get("GO domain", "")
        evidence = row.get("GO term evidence code", "")

        if pd.isna(gene_id) or str(gene_id).strip() in ("", "nan"):
            continue

        key = normalise_go_name(go_name)

        if key not in gmt:
            gmt[key] = {
                "go_domain":  str(go_dom).strip(),
                "go_id":      str(go_id).strip(),
                "genes":      [],
                "definition": str(go_def).strip(),
                "evidence":   set(),
                "raw_name":   go_name
            }

        gene_str = str(gene_id).strip()
        if gene_str not in gmt[key]["genes"]:
            gmt[key]["genes"].append(gene_str)

        gmt[key]["evidence"].add(str(evidence).strip())

    for key in gmt:
        gmt[key]["evidence"] = sorted(gmt[key]["evidence"])

    return gmt


def write_gmt(gmt: dict, out_path: str, label: str):
    """
    Format per row:
    GO_TERM_NAME <tab> GO_DOMAIN <tab> gene1 <tab> gene2 ...

    Column 2 is GO domain (biological_process, molecular_function, cellular_component)
    """
    written = skipped = 0
    with open(out_path, "w") as fh:
        for go_key, data in sorted(gmt.items()):
            if not data["genes"]:
                skipped += 1
                continue
            parts = [go_key, data["go_domain"]] + data["genes"]
            fh.write("\t".join(parts) + "\n")
            written += 1

    print(f"✓ {label}")
    print(f"  -> {out_path}")
    print(f"  {written} gene sets | {skipped} empty sets skipped")


def write_json(gmt: dict, out_path: str):
    records = []
    for go_key, data in sorted(gmt.items()):
        records.append({
            "gene_set_name":      go_key,
            "raw_go_name":        data["raw_name"],
            "go_id":              data["go_id"],
            "go_domain":          data["go_domain"],
            "evidence_codes":     data["evidence"],
            "go_term_definition": data["definition"],
            "gene_count":         len(data["genes"])
        })

    with open(out_path, "w") as fh:
        json.dump(records, fh, indent=2)

    print(f"✓ JSON written -> {out_path}")
    print(f"  {len(records)} GO term entries")


def print_preview(gmt: dict, label: str, n: int = 3):
    print(f"\n-- {label} preview (first {n} sets) --")
    for i, (key, data) in enumerate(sorted(gmt.items())):
        if i >= n:
            break
        genes = data["genes"][:3]
        suffix = f" ... ({len(data['genes'])} total)" if len(data["genes"]) > 3 else ""
        print(f"  {key}\t{data['go_domain']}\t{chr(9).join(genes)}{suffix}")
    print()




def main():
    path = sys.argv[1] if len(sys.argv) > 1 else INPUT_FILE

    if not Path(path).exists():
        print(f"File not found: {path}")
        print(f"Usage: python biomart_to_gmt.py your_biomart_file.tsv")
        sys.exit(1)

    print(f"\nLoading {path} ...")
    df = load_and_validate(path)
    print(f"  {len(df):,} rows | "
          f"{df['Gene stable ID'].nunique():,} unique genes | "
          f"{df['GO term name'].nunique():,} unique GO terms\n")

    # ── build all three GMT dicts ─────────────────────────────────────────────
    print("Building gene sets...")
    gmt_gene_stable = build_gmt(df, "Gene stable ID")
    gmt_transcript  = build_gmt(df, "Transcript stable ID")
    gmt_gene_name   = build_gmt(df, "Gene name")

    # ── previews ──────────────────────────────────────────────────────────────
    print_preview(gmt_gene_stable, "Gene stable ID GMT")
    print_preview(gmt_transcript,  "Transcript stable ID GMT")
    print_preview(gmt_gene_name,   "Gene name GMT")

    # ── write outputs ─────────────────────────────────────────────────────────
    write_gmt(gmt_gene_stable, GMT_GENE_STABLE, "Gene stable ID GMT")
    print()
    write_gmt(gmt_transcript,  GMT_TRANSCRIPT,  "Transcript stable ID GMT")
    print()
    write_gmt(gmt_gene_name,   GMT_GENE_NAME,   "Gene name GMT")
    print()
    write_json(gmt_gene_stable, JSON_OUT)

    print(f"\nDone. Output files:")
    print(f"  {GMT_GENE_STABLE}")
    print(f"  {GMT_TRANSCRIPT}")
    print(f"  {GMT_GENE_NAME}")
    print(f"  {JSON_OUT}")


if __name__ == "__main__":
    main()