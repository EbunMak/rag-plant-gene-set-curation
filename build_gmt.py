"""
build_gmt.py
============
Builds two GMT files from the gene conversion pipeline:

  1. reconstructed.gmt   — gene sets rebuilt from converted + validated genes
  2. original_subset.gmt — matching rows from the original GMT, renamed to
                           match the reconstructed names

Name normalisation for matching (applied to both sides):
  lowercase + strip all '-' and '_'

Checks directory
  One JSON file per TraesCS ID, e.g. TraesCS7B02G187400.json
  { "Gene": "...", "Validation": "yes"|"no", "Supporting Extract": "...", "PMIDS": [...] }
  The filename stem is the authoritative TraesCS ID (the Gene field may contain typos).
  Only genes with Validation == "yes" are included.
  If the checks directory is absent or empty, all converted entries are used.

Usage
  python build_gmt.py \
      --results  conversion_results.json \
      --checks   checks/ \
      --original original.gmt \
      --out-reconstructed reconstructed.gmt \
      --out-subset        original_subset.gmt
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

CONVERTED_TAGS = {"converted", "converted_normalised"}


def normalize(name: str) -> str:
    """Lowercase and strip all hyphens and underscores."""
    return re.sub(r"[-_]", " ", name.lower())


def load_checks(checks_path: Path) -> dict[str, set[str]]:
    """
    Returns { gene_set_name: {traescs_id, ...} } for validated genes only.

    Expected structure:
        checks_dir/
          <model>/                          e.g. deepseek-r1:8b/
            <gene_set_name>/               e.g. 'de novo' cotranslational protein folding/
              TraesCS1A02G420500.json
              TraesCS7B02G187400.json
              ...
    """
    if not checks_path.exists():
        print(f"[warn] checks path not found: {checks_path} — skipping validation filter")
        return {}

    validated: dict[str, set[str]] = defaultdict(set)
    total = skipped = 0

    for f in checks_path.rglob("*.json"):
        # Must be exactly 3 levels deep: model / gene_set / TraesCS.json
        try:
            rel_parts = f.relative_to(checks_path).parts
        except ValueError:
            continue
        # if len(rel_parts) != 3:
        #     continue

        gene_set, _ = rel_parts
        traescs_id = f.stem   # filename is authoritative — avoids Gene field typos

        try:
            data = json.loads(f.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"[warn] skipping {f}: {e}")
            skipped += 1
            continue

        if str(data.get("Validation", "")).strip().lower() == "yes":
            validated[gene_set].add(traescs_id)
            total += 1

    print(
        f"[info] {total} validated gene-geneset pairs "
        f"across {len(validated)} gene sets "
        f"({skipped} files skipped)"
    )
    return dict(validated)

def load_all_checked(checks_path: Path) -> dict[str, set[str]]:
    """
    Returns { gene_set_name: {traescs_id, ...} } for ALL checked genes,
    regardless of Validation value. Used for the --require-full-coverage filter.
    """
    if not checks_path.exists():
        return {}

    all_checked: dict[str, set[str]] = defaultdict(set)

    for f in checks_path.rglob("*.json"):
        try:
            rel_parts = f.relative_to(checks_path).parts
        except ValueError:
            continue
        gene_set = rel_parts[0]
        all_checked[gene_set].add(f.stem)

    return dict(all_checked)


def load_results(results_path: Path) -> list[dict]:
    data = json.loads(results_path.read_text(encoding="utf-8"))
    print(f"[info] {len(data)} entries in conversion_results.json")
    return data


def build_reconstructed(
    results: list[dict],
    validated: dict[str, set[str]],   # { gene_set: {traescs_id} }, empty = no filter
) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = defaultdict(set)
    skipped_tag = 0

    # Index all converted genes by gene set
    converted_by_geneset: dict[str, set[str]] = defaultdict(set)
    for r in results:
        if r.get("tag") not in CONVERTED_TAGS:
            skipped_tag += 1
            continue
        traescs_id = (r.get("traescs_id") or "").strip()
        gene_set   = (r.get("gene_set")   or "").strip()
        if traescs_id and gene_set:
            converted_by_geneset[gene_set].add(traescs_id.split(".")[0])

    if validated:
        # Union of gene sets seen in either source
        all_gene_sets = set(converted_by_geneset) | set(validated)
        skipped_no_genes = 0

        for gs in all_gene_sets:
            val_ids   = validated.get(gs, set())
            converted = converted_by_geneset.get(gs, set())

            if val_ids:
                # Validated IDs are authoritative — use them as the gene list.
                # Converted-but-not-validated genes are excluded when a checks
                # directory is present.
                gene_sets[gs] = val_ids
            elif converted:
                # Gene set was generated and converted but not yet checked —
                # include converted genes so the gene set isn't silently dropped.
                gene_sets[gs] = converted
            else:
                skipped_no_genes += 1

        print(
            f"[info] {len(gene_sets)} reconstructed gene sets | "
            f"skipped {skipped_tag} non-converted | "
            f"{skipped_no_genes} gene sets had neither converted nor validated genes"
        )
    else:
        # No validation filter — use all converted gene sets as-is
        gene_sets = converted_by_geneset
        print(
            f"[info] {len(gene_sets)} reconstructed gene sets | "
            f"skipped {skipped_tag} non-converted"
        )

    return dict(gene_sets)

def write_reconstructed_gmt(
    gene_sets: dict[str, set[str]],
    out_path: Path,
    descriptions: dict[str, str] | None = None,
):
    lines = 0
    with out_path.open("w", encoding="utf-8") as f:
        for name in sorted(gene_sets):
            genes = sorted(gene_sets[name])
            desc = (descriptions or {}).get(name, "")
            f.write(f"{name}\t{desc}\t" + "\t".join(genes) + "\n")
            lines += 1
    print(f"[info] reconstructed GMT → {out_path} ({lines} gene sets)")


def build_original_subset(
    original_gmt: Path,
    reconstructed_names: dict[str, str],  # norm → reconstructed name
) -> tuple[list[str], dict[str, str]]:
    """
    Read original GMT; for each row whose normalised name matches a
    reconstructed gene set, replace the name with the reconstructed name.
    Returns (GMT lines, {reconstructed_name: description_column}).
    """
    if not original_gmt.exists():
        print(f"[warn] original GMT not found: {original_gmt}")
        return [], {}

    subset: list[str] = []
    descriptions: dict[str, str] = {}
    matched = unmatched = 0

    with original_gmt.open(encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            orig_name = parts[0]
            norm = normalize(orig_name)

            if norm in reconstructed_names:
                recon_name = reconstructed_names[norm]
                parts[0] = recon_name
                subset.append("\t".join(parts))
                # second column is the description (GO domain / label)
                if len(parts) > 1:
                    descriptions[recon_name] = parts[1]
                matched += 1
            else:
                unmatched += 1

    print(
        f"[info] original GMT: {matched} matched / "
        f"{unmatched} unmatched gene sets"
    )
    return subset, descriptions


def write_subset_gmt(lines: list[str], out_path: Path):
    with out_path.open("w", encoding="utf-8") as f:
        for line in lines:
            f.write(line + "\n")
    print(f"[info] original subset GMT → {out_path} ({len(lines)} gene sets)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Build GMTs from conversion pipeline output")
    parser.add_argument("--results",          default="conversion_results.json",
                        help="Path to conversion_results.json")
    parser.add_argument("--checks",           default="checks/",
                        help="Path to checks directory (or single JSON file)")
    parser.add_argument("--original",         default="original.gmt",
                        help="Path to original GMT file")
    parser.add_argument("--out-reconstructed", default="out/reconstructed.gmt",
                        help="Output path for reconstructed GMT")
    parser.add_argument("--out-subset",       default="out/original_subset.gmt",
                        help="Output path for original GMT subset")
    parser.add_argument("--require-full-coverage", action="store_true",
                        help="Only include gene sets where every converted gene "
                             "has a check file (coverage = 100%%), regardless of "
                             "validation result")
    args = parser.parse_args()

    results_path     = Path(args.results)
    checks_path      = Path(args.checks)
    original_gmt     = Path(args.original)
    out_reconstructed = Path(args.out_reconstructed)
    out_subset        = Path(args.out_subset)

    if not results_path.exists():
        sys.exit(f"[error] results file not found: {results_path}")

    # 1. Load checks
    validated_ids = load_checks(checks_path)

    # 2. Load conversion results
    results = load_results(results_path)

    # 3. Build reconstructed gene sets
    gene_sets = build_reconstructed(results, validated_ids)

    if not gene_sets:
        sys.exit("[error] no gene sets produced — check tags and validation")

    # 3b. Optionally filter to gene sets with 100% check coverage
    if args.require_full_coverage:
        all_checked = load_all_checked(checks_path)
        before = len(gene_sets)
        fully_covered = {}
        for gs, genes in gene_sets.items():
            checked_ids = all_checked.get(gs, set())
            unchecked = genes - checked_ids
            if not unchecked:
                fully_covered[gs] = genes
        gene_sets = fully_covered
        dropped_cov = before - len(gene_sets)
        print(
            f"[info] --require-full-coverage: kept {len(gene_sets)} / {before} gene sets "
            f"({dropped_cov} dropped due to incomplete check coverage)"
        )

    # 4. Build normalised name → reconstructed name lookup
    reconstructed_norm = {normalize(name): name for name in gene_sets}

    # 5. Extract matching rows from original GMT (determines the valid name set)
    subset_lines, descriptions = build_original_subset(original_gmt, reconstructed_norm)

    # 6. Filter gene_sets to only those matched in the original GMT
    matched_names = {line.split("\t")[0] for line in subset_lines}
    gene_sets_matched = {name: genes for name, genes in gene_sets.items()
                         if name in matched_names}

    dropped = len(gene_sets) - len(gene_sets_matched)
    if dropped:
        print(f"[info] dropped {dropped} reconstructed gene sets not found in original GMT")

    # 7. Write reconstructed GMT (matched only)
    write_reconstructed_gmt(gene_sets_matched, out_reconstructed, descriptions)

    # 8. Write original subset GMT
    if subset_lines:
        write_subset_gmt(subset_lines, out_subset)
    else:
        print("[warn] no matching gene sets found between original GMT and reconstructed sets")

    # 9. Summary
    print("\n--- Summary ---")
    print(f"  Reconstructed gene sets : {len(gene_sets_matched)}")
    print(f"  Original subset matched : {len(subset_lines)}")
    total_genes = sum(len(v) for v in gene_sets_matched.values())
    print(f"  Total TraesCS IDs       : {total_genes}")


if __name__ == "__main__":
    main()