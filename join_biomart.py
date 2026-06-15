"""
merge_biomart.py
----------------
Memory-efficient merge of all BioMart export files.
Builds the master table using a dict keyed by Gene stable ID
instead of pandas outer merges, avoiding RAM explosion.
"""

import os
import csv
import json
from collections import defaultdict

BIOMART_DIR = "data/mart_exports"
OUTPUT_FILE = "data/merged_biomart.tsv"


def detect_separator(file_path: str) -> str:
    with open(file_path, "r", encoding="utf-8") as f:
        first_line = f.readline()
    return "\t" if "\t" in first_line else ","


def merge_biomart_files():
    all_files = sorted([
        f for f in os.listdir(BIOMART_DIR)
        if f.endswith(".txt") or f.endswith(".tsv")
    ])

    if not all_files:
        print(f"No .txt or .tsv files found in {BIOMART_DIR}")
        return

    print(f"Found {len(all_files)} files\n")

    # master[gene_id] = {col: value, col: value, ...}
    master = defaultdict(dict)
    all_columns = []   # preserves column order across files

    for filename in all_files:
        if filename == "mart_export.tsv":
            print(f"  Skipping existing {filename}")
            continue
        file_path = os.path.join(BIOMART_DIR, filename)
        sep = detect_separator(file_path)

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter=sep)
                headers = reader.fieldnames

                if not headers or "Gene stable ID" not in headers:
                    print(f"  Skipping {filename} — no 'Gene stable ID' column")
                    continue

                # register new columns in order
                for col in headers:
                    if col != "Gene stable ID" and col not in all_columns:
                        all_columns.append(col)

                rows = 0
                for row in reader:
                    gene_id = row.get("Gene stable ID", "").strip()
                    if not gene_id or gene_id.startswith("ENSRNA"):
                        continue   # skip non-coding RNA entries
                    for col in headers:
                        if col == "Gene stable ID":
                            continue
                        val = row.get(col, "").strip()
                        # keep existing value if already populated
                        if val and col not in master[gene_id]:
                            master[gene_id][col] = val
                    rows += 1

            print(f"  {filename}: {rows:,} rows, "
                  f"new cols: {[c for c in headers if c != 'Gene stable ID']}")

        except Exception as e:
            print(f"  ERROR in {filename}: {e}")
            continue

    if not master:
        print("No data loaded.")
        return

    print(f"\nWriting {len(master):,} genes, {len(all_columns)} columns...")

    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    fieldnames = ["Gene stable ID"] + all_columns

    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for gene_id in sorted(master.keys()):
            row = {"Gene stable ID": gene_id}
            row.update(master[gene_id])
            writer.writerow(row)

    print(f"Saved to {OUTPUT_FILE}")
    print(f"Columns: {fieldnames}")


if __name__ == "__main__":
    merge_biomart_files()