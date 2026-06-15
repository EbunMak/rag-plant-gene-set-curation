
import pandas as pd
import json
import os
from biomart_to_gmt import load_and_validate

BIOMART_EXPORT_TSV = "data/mart_exports/mart_export.tsv"
# convert gene stable IDs to gene names and transcript stable IDs using TSV from BioMart export
def build_id_lookup(df: pd.DataFrame) -> dict:
    """
    Build a lookup dict keyed by Transcript stable ID.
    Call once at load time, then query instantly for any transcript.

    Usage:
        lookup = build_id_lookup(df)
        info   = lookup.get("TraesCS7D02G315900.1")
        print(info["gene_stable_id"])   # TraesCS7D02G315900
        print(info["gene_name"])        # WM12
        print(info["gene_synonym"])     # LOC100136962
    """
    lookup = {}

    for _, row in df.iterrows():
        transcript_id = str(row["Transcript stable ID"]).strip()
        if transcript_id in lookup:
            continue                    # already stored from first occurrence
        lookup[transcript_id] = {
            "gene_stable_id":   str(row["Gene stable ID"]).strip(),
            "gene_name":        str(row["Gene name"]).strip(),
            "transcript_name":  str(row["Transcript stable ID"]).strip()
        }

    return lookup


def convert_gene_iD(gene_id: str) -> dict:
    """
    Convert a gene ID (stable or transcript) to all related IDs using the lookup.
    Returns a dict with keys: gene_stable_id, gene_name, transcript_name, gene_synonym.
    If not found, returns None.
    """
    df = load_and_validate(BIOMART_EXPORT_TSV)
    lookup = build_id_lookup(df)
    # First try as transcript ID
    if gene_id in lookup:
        return lookup[gene_id]

    # Not found
    return None