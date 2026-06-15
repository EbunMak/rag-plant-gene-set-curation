"""
subset_direct_gmt.py
--------------------
Filters a direct-prompting GMT to only the gene sets that appear in the
reconstructed GMT (matched by normalised name).

Usage:
    python subset_direct_gmt.py \
        --reconstructed out/reconstructed.gmt \
        --direct        out/conversion_direct/deepseek-r1:8b/direct_prompt_curated.gmt \
        --out           out/conversion_direct/deepseek-r1:8b/direct_prompt_subset.gmt
"""

import argparse
import re
from pathlib import Path


def normalize(name: str) -> str:
    return re.sub(r"[-_]", " ", name.lower())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reconstructed", required=True,
                        help="Reconstructed GMT (the 732-matched one)")
    parser.add_argument("--direct",        required=True,
                        help="Direct-prompting GMT to subset")
    parser.add_argument("--out",           required=True,
                        help="Output path for the subsetted direct GMT")
    args = parser.parse_args()

    reconstructed = Path(args.reconstructed)
    direct        = Path(args.direct)
    out_path      = Path(args.out)

    # Load normalised names from reconstructed GMT
    recon_norms: set[str] = set()
    with reconstructed.open(encoding="utf-8") as f:
        for line in f:
            name = line.split("\t")[0].strip()
            if name:
                recon_norms.add(normalize(name))
    print(f"Reconstructed GMT : {len(recon_norms)} gene sets")

    # Filter direct GMT
    matched = unmatched = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with direct.open(encoding="utf-8") as fin, out_path.open("w", encoding="utf-8") as fout:
        for line in fin:
            name = line.split("\t")[0].strip()
            if not name:
                continue
            if normalize(name) in recon_norms:
                fout.write(line if line.endswith("\n") else line + "\n")
                matched += 1
            else:
                unmatched += 1

    print(f"Direct GMT        : {matched} matched, {unmatched} not in reconstructed")
    print(f"Output            : {out_path}")


if __name__ == "__main__":
    main()
