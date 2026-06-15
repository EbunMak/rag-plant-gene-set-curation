"""
match_gmt_versions.py
=====================
Matches gene sets between a v1.1 GMT and a v2.1 GMT by normalised name,
replacing v1.1 gene lists with v2.1 gene lists where a match is found.

Name normalisation: lowercase + strip all '-' and '_'

Outputs:
  - matched_v21.gmt        : gene sets from v1.1 that matched v2.1, using v2.1 genes
  - unmatched_v11.txt      : gene sets in v1.1 with no v2.1 equivalent
  - new_in_v21.txt         : gene sets in v2.1 not present in v1.1

Usage:
    python match_gmt_versions.py --v11 original_v11.gmt --v21 original_v21.gmt
"""

import argparse
import re
from pathlib import Path


def normalise(name: str) -> str:
    return re.sub(r"[-_]", "", name.lower())


def load_gmt(path: Path) -> dict[str, tuple[str, list[str]]]:
    """
    Returns { raw_name: (description, [genes]) }
    """
    gmt = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0].strip():
                continue
            name  = parts[0].strip()
            desc  = parts[1].strip() if len(parts) > 1 else name
            genes = [g.strip() for g in parts[2:] if g.strip()]
            gmt[name] = (desc, genes)
    return gmt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v11", required=True, help="v1.1 GMT file")
    parser.add_argument("--v21", required=True, help="v2.1 GMT file")
    parser.add_argument("--out", default="matched_v21.gmt")
    args = parser.parse_args()

    v11_path = Path(args.v11)
    v21_path = Path(args.v21)
    out_path = Path(args.out)

    v11 = load_gmt(v11_path)
    v21 = load_gmt(v21_path)

    print(f"v1.1 GMT : {len(v11):,} gene sets")
    print(f"v2.1 GMT : {len(v21):,} gene sets")

    # Build normalised name → v2.1 raw name lookup
    v21_norm = {normalise(name): name for name in v21}

    matched   = []   # (v11_name, v21_name, v21_genes)
    unmatched = []   # v11 names with no v2.1 equivalent

    for v11_name, (desc, _) in v11.items():
        norm = normalise(v11_name)
        if norm in v21_norm:
            v21_name  = v21_norm[norm]
            v21_desc, v21_genes = v21[v21_name]
            matched.append((v11_name, v21_name, v21_genes))
        else:
            unmatched.append(v11_name)

    # Gene sets in v2.1 not in v1.1
    v11_norms = {normalise(n) for n in v11}
    new_in_v21 = [name for norm, name in v21_norm.items() if norm not in v11_norms]

    # ── Write matched GMT ────────────────────────────────────────────────────
    with open(out_path, "w") as f:
        for v11_name, v21_name, genes in matched:
            f.write(f"{v21_name}\t{v21_name}\t" + "\t".join(genes) + "\n")

    # ── Write unmatched v1.1 names ───────────────────────────────────────────
    unmatched_path = out_path.with_stem(out_path.stem + "_unmatched_v11")
    unmatched_path = unmatched_path.with_suffix(".txt")
    with open(unmatched_path, "w") as f:
        f.write("\n".join(unmatched) + "\n")

    # ── Write new v2.1 names not in v1.1 ────────────────────────────────────
    new_path = out_path.with_stem(out_path.stem + "_new_in_v21")
    new_path = new_path.with_suffix(".txt")
    with open(new_path, "w") as f:
        f.write("\n".join(new_in_v21) + "\n")

    # ── Summary ──────────────────────────────────────────────────────────────
    print(f"\nMatched        : {len(matched):,}")
    print(f"Unmatched v1.1 : {len(unmatched):,}  (in v1.1 but not v2.1)")
    print(f"New in v2.1    : {len(new_in_v21):,}  (in v2.1 but not v1.1)")
    print(f"\nOutputs:")
    print(f"  {out_path}")
    print(f"  {unmatched_path}")
    print(f"  {new_path}")


if __name__ == "__main__":
    main()