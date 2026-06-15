#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np

def analyze_gene_sets(csv_path, output_dir):
    """Analyze gene set comparison CSV and save statistics + gene lists"""
    
    # Load comparison table
    df = pd.read_csv(csv_path)
    
    # Ensure numerical columns exist and are of correct type
    for col in ["# Common", "# New", "# Original"]:
        if col not in df.columns:
            df[col] = 0
    df[["# Common", "# New", "# Original"]] = df[["# Common", "# New", "# Original"]].fillna(0).astype(int)
    
    # Compute # of genes in original but missing in new
    df["# Lost"] = df["# Original"] - df["# Common"]
    
    # Capture ALL print output
    stats_output = []
    
    def print_stats(msg):
        print(msg)
        stats_output.append(msg)
    
    print_stats("\n=== BASIC DATA ===")
    print_stats(f"Total gene sets compared: {len(df)}")
    
    print_stats("\n=== OVERLAP STATISTICS (# COMMON GENES) ===")
    print_stats(str(df["# Common"].describe()))
    
    print_stats("\nGene sets with NO common genes:")
    no_common = df[df["# Common"] == 0]["Gene Set Name"].tolist()
    print_stats(f"  Count: {len(no_common)}")
    print_stats(f"  Names: {no_common}")
    
    print_stats("\nGene sets with HIGH overlap (top 10):")
    top_overlap = df.sort_values("# Common", ascending=False)[["Gene Set Name", "# Common"]].head(10)
    print_stats(str(top_overlap))
    
    print_stats("\n=== NEW GENE STATISTICS (# NEW GENES) ===")
    print_stats(str(df["# New"].describe()))
    
    print_stats("\nGene sets with NO new genes added (fully consistent):")
    no_new = df[df["# New"] == 0]["Gene Set Name"].tolist()
    print_stats(f"  Count: {len(no_new)}")
    print_stats(f"  Names: {no_new}")
    
    print_stats("\nGene sets with MANY new genes added (top 10):")
    top_new = df.sort_values("# New", ascending=False)[["Gene Set Name", "# New"]].head(10)
    print_stats(str(top_new))
    
    print_stats("\n=== LOSS STATISTICS (# ORIGINAL GENES NOT IN NEW) ===")
    print_stats(str(df["# Lost"].describe()))
    
    print_stats("\nGene sets that lost ALL original genes:")
    total_loss = df[df["# Lost"] == df["# Original"]]["Gene Set Name"].tolist()
    print_stats(f"  Count: {len(total_loss)}")
    print_stats(f"  Names: {total_loss}")
    
    print_stats("\nGene sets that retained MOST original genes (lowest loss, top 10):")
    low_loss = df.sort_values("# Lost", ascending=True)[["Gene Set Name", "# Lost"]].head(10)
    print_stats(str(low_loss))
    
    # GLOBAL DATABASE-LEVEL STATISTICS
    def parse_gene_list(col):
        if isinstance(col, str) and col.strip():
            return set(g.strip() for g in col.split(",") if g.strip())
        return set()
    
    # Parse existing columns
    df["Common Set"] = df["Common Genes"].apply(parse_gene_list)
    df["Lost Set"] = df["Lost Genes"].apply(parse_gene_list)
    
    if "Original Genes" in df.columns:
        df["Original Set"] = df["Original Genes"].apply(parse_gene_list)
    else:
        df["Original Set"] = df.apply(lambda row: row["Common Set"].union(row["Lost Set"]), axis=1)
    
    df["New Set"] = df.apply(
        lambda row: row["Common Set"].union(parse_gene_list(row["Newly Added Genes"])),
        axis=1
    )
    
    # Global database unions
    global_original = set().union(*df["Original Set"].tolist())
    global_new = set().union(*df["New Set"].tolist())
    global_lost = global_original - global_new
    global_gained = global_new - global_original
    
    print_stats("\n=== GLOBAL DATABASE-LEVEL GENE DIFFERENCES ===")
    print_stats(f"Total unique genes in ORIGINAL DB: {len(global_original)}")
    print_stats(f"Total unique genes in NEW DB: {len(global_new)}")
    print_stats(f"Genes LOST from original DB: {len(global_lost)}")
    print_stats(f"Genes GAINED in new DB: {len(global_gained)}")
    
    # SAVE FILES IN SAME DIR AS INPUT CSV
    os.makedirs(output_dir, exist_ok=True)
    
    # Save statistics to txt
    stats_file = os.path.join(output_dir, "statistics.txt")
    with open(stats_file, "w") as f:
        f.write("\n".join(stats_output))
    
    # Save gene lists
    pd.DataFrame({"Genes Lost": list(global_lost)}).to_csv(
        os.path.join(output_dir, "genes_lost_globally.csv"), index=False)
    pd.DataFrame({"Genes Gained": list(global_gained)}).to_csv(
        os.path.join(output_dir, "genes_gained_globally.csv"), index=False)
    
    print(f"\n📁 Saved to {output_dir}:")
    print(f"  • statistics.txt")
    print(f"  • genes_lost_globally.csv") 
    print(f"  • genes_gained_globally.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze gene set comparison CSV")
    parser.add_argument("csv_path", help="Path to gene_set_comparison.csv")
    parser.add_argument("--config", required=True, help="Configuration name (for output dir)")
    
    args = parser.parse_args()
    
    # Create output dir in same location as input CSV
    input_dir = os.path.dirname(os.path.abspath(args.csv_path))
    output_dir = os.path.join(input_dir, args.config)
    
    analyze_gene_sets(args.csv_path, output_dir)
