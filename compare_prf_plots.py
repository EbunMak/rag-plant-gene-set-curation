import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


def add_summary_stats(df):
    metrics = ["Precision", "Recall", "F1"]
    for metric in metrics:
        print("\n=== {} ===".format(metric))
        stats = df.groupby("Configuration")[metric].agg(["mean", "median"])
        print(stats.round(3))


def get_configuration_order(df):
    f1_medians = df.groupby("Configuration")["F1"].median().sort_values(ascending=False)
    return list(f1_medians.index)


def get_palette(config_order):
    base_palette = sns.color_palette("Set2", n_colors=len(config_order))
    return {cfg: color for cfg, color in zip(config_order, base_palette)}


def plot_single_metric(ax, df, metric, config_order, palette_map, panel_label):
    """Plot one metric panel in the narrow combined style."""

    # Grey background panel
    ax.set_facecolor("#ebebeb")

    sns.boxplot(
        x="Configuration",
        y=metric,
        data=df,
        order=config_order,
        palette=[palette_map[cfg] for cfg in config_order],
        showfliers=False,
        notch=False,
        ax=ax,
        linewidth=1.2,
        boxprops=dict(alpha=0.7),
        width=0.45,
    )

    # Median and mean labels
    grouped = df.groupby("Configuration")[metric]
    medians = grouped.median()
    means   = grouped.mean()
    for i, cfg in enumerate(config_order):
        ax.text(i, 1.07, f"med={medians[cfg]:.2f}",
                ha="center", va="bottom", fontsize=7,
                fontweight="bold", color="#111111")
        ax.text(i, 1.02, f"mean={means[cfg]:.2f}",
                ha="center", va="bottom", fontsize=6.5,
                fontweight="normal", color="#555555")

    ax.set_ylim(-0.05, 1.20)
    ax.set_ylabel(metric, fontsize=10, fontweight="bold", labelpad=8)
    ax.set_xlabel("")

    ax.set_xticks(range(len(config_order)))
    ax.set_xticklabels(config_order, rotation=35, ha="right",
                       fontsize=8, fontweight="bold")
    ax.tick_params(axis="y", labelsize=8)
    for label in ax.get_yticklabels():
        label.set_fontweight("bold")

    # Light grid lines like the reference plot
    ax.grid(True, axis="y", linestyle="-", linewidth=0.5, alpha=0.5, color="white")
    ax.set_axisbelow(True)

    # Panel label (a, b, c)
    ax.text(-0.22, 1.18, panel_label, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top")


def plot_combined(df, out_dir, config_order, palette_map):
    """Single narrow figure with three side-by-side panels (a, b, c)."""
    os.makedirs(out_dir, exist_ok=True)

    plt.rcParams["font.family"] = "Liberation Serif"
    plt.rcParams["font.weight"] = "bold"

    metrics = ["Precision", "Recall", "F1"]
    labels  = ["a)", "b)", "c)"]

    # Narrower figure + tighter spacing between panels
    fig, axes = plt.subplots(1, 3, figsize=(7, 3.8), sharey=False)

    for ax, metric, label in zip(axes, metrics, labels):
        plot_single_metric(ax, df, metric, config_order, palette_map, label)

    fig.text(0.5, -0.04, "Configuration", ha="center",
             fontsize=10, fontweight="bold")

    plt.subplots_adjust(wspace=0.45)
    plt.tight_layout(rect=[0, 0, 1, 1])
    out_path = os.path.join(out_dir, "prf_combined.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved combined PRF plot to {out_path}")


def load_prf_tables(csv_paths, labels):
    all_dfs = []
    for path, label in zip(csv_paths, labels):
        df = pd.read_csv(path)
        for col in ["Precision", "Recall", "F1"]:
            if col not in df.columns:
                raise ValueError("{} not found in {}".format(col, path))
        df["Configuration"] = label
        all_dfs.append(df)
    return pd.concat(all_dfs, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(
        description="Compare multiple PRF tables via boxplots."
    )
    parser.add_argument("--csvs", nargs="+", required=True)
    parser.add_argument("--labels", nargs="+", required=True)
    parser.add_argument("--out_dir", type=str, default="prf_comparison_plots")
    args = parser.parse_args()

    if len(args.csvs) != len(args.labels):
        raise ValueError("Number of --csvs must match number of --labels")

    df = load_prf_tables(args.csvs, args.labels)
    add_summary_stats(df)

    config_order = get_configuration_order(df)
    palette_map = get_palette(config_order)

    plot_combined(df, args.out_dir, config_order, palette_map)


if __name__ == "__main__":
    main()