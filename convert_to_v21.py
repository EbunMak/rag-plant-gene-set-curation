"""
convert_to_v21.py
=================
Downloads the IWGSC all-correspondances file from URGI and uses it to
build a v1.1 → v2.1 TraesCS ID mapping, then applies it to:

  1. conversion_results.json  → adds 'traescs_v21_id' field per entry
  2. reconstructed.gmt        → writes reconstructed_v21.gmt
  3. original_subset.gmt      → writes original_subset_v21.gmt

File format (space-separated, 4 columns):
  v1.0                  v1.1                  v2.1                  css2014
  TraesCS1A01G000100LC  TraesCS1A02G000100LC  TraesCS1A03G0000100LC -
  TraesCS1A01G000100    TraesCS1A02G000100    TraesCS1A03G0000200   -

Usage:
    python convert_to_v21.py
    python convert_to_v21.py --results conversion_results.json \
                             --gmt-reconstructed reconstructed.gmt \
                             --gmt-subset original_subset.gmt
"""

import argparse
import json
import re
import zipfile
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

CORRESPONDANCES_URL = (
    "https://urgi.versailles.inra.fr/download/iwgsc/"
    "IWGSC_RefSeq_Annotations/v2.1/iwgsc_refseq_all_correspondances.zip"
)
CORRESPONDANCES_ZIP = Path("iwgsc_refseq_all_correspondances.zip")


# ---------------------------------------------------------------------------
# Step 1 — Download
# ---------------------------------------------------------------------------

def download_correspondances(zip_path: Path = CORRESPONDANCES_ZIP):
    if zip_path.exists():
        print(f"[download] Already exists: {zip_path} — skipping download")
        return
    print(f"[download] Fetching {CORRESPONDANCES_URL} ...")
    r = requests.get(CORRESPONDANCES_URL, stream=True, timeout=60)
    r.raise_for_status()
    with open(zip_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=65536):
            f.write(chunk)
    print(f"[download] Saved → {zip_path} ({zip_path.stat().st_size // 1024} KB)")


# ---------------------------------------------------------------------------
# Step 2 — Parse correspondances file
# ---------------------------------------------------------------------------

def build_mapping(zip_path: Path = CORRESPONDANCES_ZIP) -> tuple[dict, dict]:
    """
    Parse the correspondances zip.

    File format (space-separated):
        v1.0                  v1.1                  v2.1                  css2014
        TraesCS1A01G000100LC  TraesCS1A02G000100LC  TraesCS1A03G0000100LC -

    Returns:
        v11_to_v21: { TraesCS...02G... : TraesCS...03G... }
        v10_to_v21: { TraesCS...01G... : TraesCS...03G... }
    """
    v11_to_v21: dict[str, str] = {}
    v10_to_v21: dict[str, str] = {}
    total = mapped_v11 = mapped_v10 = skipped = 0

    with zipfile.ZipFile(zip_path) as zf:
        data_files = [n for n in zf.namelist() if not n.endswith("/")]
        print(f"[mapping] Files in zip: {data_files}")

        for fname in data_files:
            with zf.open(fname) as raw:
                lines = raw.read().decode("utf-8", errors="replace").splitlines()

            if not lines:
                continue

            print(f"[mapping] Parsing {fname} — {len(lines):,} lines")
            print(f"[mapping] Header : {lines[0]}")
            print(f"[mapping] Row 1  : {lines[1] if len(lines) > 1 else 'N/A'}")

            for line in lines[1:]:
                line = line.strip()
                if not line:
                    continue
                total += 1

                # Columns are space-separated; split by whitespace
                # Format: v1.0_id  v1.1_id  v2.1_id  css2014(ignored)
                tokens = line.split()

                # Match by position (columns are fixed) with pattern fallback
                v10_id = tokens[0].split(".")[0] if len(tokens) > 0 else ""
                v11_id = tokens[1].split(".")[0] if len(tokens) > 1 else ""
                v21_id = tokens[2].split(".")[0] if len(tokens) > 2 else ""

                # Validate by pattern — skip if a token doesn't look like a TraesCS ID
                if not re.match(r"TraesCS\w+01G\d+", v10_id, re.I): v10_id = ""
                if not re.match(r"TraesCS\w+02G\d+", v11_id, re.I): v11_id = ""
                if not re.match(r"TraesCS\w+03G\d+", v21_id, re.I): v21_id = ""

                if v11_id and v21_id:
                    v11_to_v21[v11_id] = v21_id
                    mapped_v11 += 1
                if v10_id and v21_id:
                    v10_to_v21[v10_id] = v21_id
                    mapped_v10 += 1
                if not v21_id:
                    skipped += 1

    print(
        f"[mapping] {total:,} rows | "
        f"{mapped_v11:,} v1.1→v2.1 | "
        f"{mapped_v10:,} v1.0→v2.1 | "
        f"{skipped:,} skipped (no v2.1 ID)"
    )
    return v11_to_v21, v10_to_v21


# ---------------------------------------------------------------------------
# Step 3 — Update conversion_results.json
# ---------------------------------------------------------------------------

def update_results(results_path: Path, v11_to_v21: dict, v10_to_v21: dict):
    with open(results_path) as f:
        results = json.load(f)

    converted = not_found = 0
    for entry in results:
        raw_id  = (entry.get("traescs_id") or "").split(".")[0]
        v21_id  = v11_to_v21.get(raw_id) or v10_to_v21.get(raw_id) or ""
        entry["traescs_v21_id"] = v21_id if v21_id else None
        if v21_id:
            converted += 1
        else:
            not_found += 1

    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

    print(
        f"[results] {len(results):,} entries updated — "
        f"{converted:,} have a v2.1 ID, {not_found:,} unmapped"
    )


# ---------------------------------------------------------------------------
# Step 4 — Convert GMT files
# ---------------------------------------------------------------------------

def convert_gmt(gmt_in: Path, gmt_out: Path, v11_to_v21: dict, v10_to_v21: dict):
    if not gmt_in.exists():
        print(f"[gmt] {gmt_in} not found — skipping")
        return

    total = converted = kept = 0

    with open(gmt_in) as fin, open(gmt_out, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                fout.write(line)
                continue

            name, desc, genes_in = parts[0], parts[1], parts[2:]
            genes_out = []

            for g in genes_in:
                base  = re.sub(r"LC$", "", g.strip().split(".")[0], flags=re.I)
                v21   = v11_to_v21.get(base) or v10_to_v21.get(base)
                total += 1
                if v21:
                    genes_out.append(v21)
                    converted += 1
                else:
                    genes_out.append(g.strip())  # keep original if no mapping
                    kept += 1

            fout.write(f"{name}\t{desc}\t" + "\t".join(genes_out) + "\n")

    print(
        f"[gmt] {gmt_in.name} → {gmt_out.name} | "
        f"{converted:,} converted, {kept:,} kept as-is (no mapping)"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Convert TraesCS v1.1 → v2.1 IDs")
    parser.add_argument("--results",           default="conversion_results.json")
    parser.add_argument("--gmt-reconstructed", default="reconstructed.gmt")
    parser.add_argument("--gmt-subset",        default="original_subset.gmt")
    args = parser.parse_args()

    # 1. Download
    download_correspondances()

    # 2. Build mapping
    v11_to_v21, v10_to_v21 = build_mapping()
    if not v11_to_v21:
        raise SystemExit("[error] Empty mapping — check the correspondances file")

    # Save mapping for reuse by other scripts
    mapping_out = Path("v11_to_v21_mapping.json")
    with open(mapping_out, "w") as f:
        json.dump(v11_to_v21, f, indent=2)
    print(f"[mapping] Saved → {mapping_out} ({len(v11_to_v21):,} entries)")

    # 3. Update conversion_results.json
    results_path = Path(args.results)
    if results_path.exists():
        update_results(results_path, v11_to_v21, v10_to_v21)
    else:
        print(f"[results] {results_path} not found — skipping")

    # 4. Convert GMTs
    for gmt_name in [args.gmt_reconstructed, args.gmt_subset]:
        gmt_in  = Path(gmt_name)
        gmt_out = gmt_in.with_stem(gmt_in.stem + "_v21")
        convert_gmt(gmt_in, gmt_out, v11_to_v21, v10_to_v21)

    print("\n--- Done ---")
    print(f"  v11→v21 mapping : {len(v11_to_v21):,} entries  →  v11_to_v21_mapping.json")
    print(f"  Results updated : {args.results}")
    print(f"  GMTs written    : *_v21.gmt")


if __name__ == "__main__":
    main()