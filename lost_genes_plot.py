import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from collections import defaultdict
import matplotlib.cm as cm

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.size": 14
})

def make_plot(model_csv_pairs):
    # Load all model data
    all_data = []
    for model_name, csv_path in model_csv_pairs.items():
        df = pd.read_csv(csv_path)
        df["model"] = model_name  # Assign model name
        df["% Loss"] = (df["# Lost"] / df["# Original"]) * 100
        all_data.append(df[["model", "% Loss"]].dropna())
    
    if not all_data:
        print("No valid data found!")
        return
    
    df_combined = pd.concat(all_data, ignore_index=True)
    
    # Get unique ordered models
    models = sorted(df_combined["model"].unique())
    df_combined["model"] = pd.Categorical(df_combined["model"], categories=models, ordered=True)
    
    # Output directory (use first CSV's dir)
    first_csv = list(model_csv_pairs.values())[0]
    base_dir = os.path.dirname(first_csv)
    plot_dir = os.path.join(base_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # HORIZONTAL ARRANGEMENT: 2 rows, auto-calculated columns
    n_models = len(models)
    ncols = min(4, (n_models + 1) // 2)  # Max 4 cols, 2 rows
    nrows = (n_models + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharex=False, sharey=False)
    if nrows * ncols == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Colors matching image (blue gradient)
    colors = cm.Blues(np.linspace(0.3, 0.9, n_models))
    
    for i, model in enumerate(models):
        ax = axes[i]
        model_data = df_combined[df_combined["model"] == model]["% Loss"]
        
        if model_data.empty:
            ax.set_visible(False)
            continue
        
        # Histogram (density=True for ridges)
        ax.hist(model_data, bins=20, density=True, alpha=0.6, color=colors[i],
                edgecolor="black", linewidth=0.5, label="Histogram")
        
        # KDE (smooth density curve)
        kde = gaussian_kde(model_data, bw_method=0.3)
        x_min, x_max = model_data.min(), model_data.max()
        x_vals = np.linspace(max(0, x_min-5), x_max+5, 600)
        y_vals = kde(x_vals)
        ax.plot(x_vals, y_vals, linewidth=2.5, color=colors[i], label="Density")
        
        # Mean vertical line
        mean_val = model_data.mean()
        ax.axvline(mean_val, color="red", linestyle="--", linewidth=2,
                   label=f"Mean = {mean_val:.2f}%")
        
        # Ridge-style labels (y="Prompt Strategy" -> model)
        ax.set_title(model.replace('-', '\n'), fontsize=14, fontweight="bold", pad=15)
        ax.set_ylabel("Density", fontsize=12, fontweight="bold")
        ax.set_xlabel("% Loss", fontsize=12, fontweight="bold")
        ax.legend(fontsize=10, frameon=True)
        # set ledgend to be bold
        ax.legend().get_frame().set_edgecolor("black")
        ax.legend().get_frame().set_linewidth(0.5)
        ax.legend().get_texts()[0].set_fontweight("bold")
        ax.legend().get_texts()[1].set_fontweight("bold")
        ax.legend().get_texts()[2].set_fontweight("bold")
        ax.grid(alpha=0.15)
        ax.tick_params(axis="both", labelsize=10)
        # set tick labels to be bold
        for label in ax.get_xticklabels():
            label.set_fontweight("bold")
        for label in ax.get_yticklabels():
            label.set_fontweight("bold")
        ax.grid(alpha=0.15)
        
        # Auto-scale each subplot independently
        ax.margins(x=0.05, y=0.1)
    
    # Hide empty subplots
    for j in range(len(models), len(axes)):
        axes[j].set_visible(False)
    
    plt.suptitle("Distribution of % Loss Across Models", fontsize=16, fontweight="bold", y=0.98)
    plt.tight_layout()
    
    # Save
    base_name = os.path.join(plot_dir, "percent_loss_horizontal")
    plt.savefig(f"{base_name}.png", dpi=400, bbox_inches="tight")
    plt.savefig(f"{base_name}.pdf", bbox_inches="tight")
    plt.savefig(f"{base_name}.svg", bbox_inches="tight")
    
    print("\nSaved horizontal plots:")
    for ext in [".png", ".pdf", ".svg"]:
        print(f" - {base_name}{ext}")
    print()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot % Loss distributions for multiple models")
    parser.add_argument("--model_csv", nargs="*", action="append", default=[],
                       help="Model=CSV pairs, e.g. '--model_csv model1=out1.csv model2=out2.csv'")
    parser.add_argument("--comparison_csv", type=str, default=None,
                       help="Fallback single CSV path (original behavior)")
    
    args = parser.parse_args()
    
    model_csv_pairs = {}
    
    # Parse model_csv pairs (e.g., model1=out1.csv model2=out2.csv)
    if args.model_csv:
        for pair_list in args.model_csv:
            for pair in pair_list:
                if "=" in pair:
                    model, csv_path = pair.split("=", 1)
                    model_csv_pairs[model.strip()] = csv_path.strip()
    
    # Fallback: single CSV (original behavior, assign default model name)
    if args.comparison_csv and os.path.exists(args.comparison_csv):
        if not model_csv_pairs:
            model_csv_pairs["Model"] = args.comparison_csv
        else:
            print(f"Warning: --comparison_csv provided but --model_csv used: {model_csv_pairs}")
    
    if not model_csv_pairs:
        print("Error: Provide --model_csv 'model1=out1.csv model2=out2.csv' or --comparison_csv")
        exit(1)
    
    make_plot(model_csv_pairs)
