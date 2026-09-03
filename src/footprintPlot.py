#!/usr/bin/env python3
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main():
    parser = argparse.ArgumentParser(
        description="Generate Differential Binding Heatmap from BINDetect results."
    )
    parser.add_argument(
        "--indir",
        type=str,
        default="bindetect_results",
        help="Input bindetect_results directory",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="summary_plots",
        help="Output directory for plots",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="bindetect",
        help="Prefix for saved plot filenames",
    )
    parser.add_argument(
        "-n",
        "--n",
        type=int,
        default=30,
        help="Number of top dynamic TFs to plot in the heatmap (default: 30)",
    )
    parser.add_argument(
        "--pvalue",
        type=float,
        default=0.05,
        help="Maximum p-value threshold for statistical significance (default: 0.05)",
    )
    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    results_file = indir / "bindetect_results.txt"

    if not results_file.exists():
        raise FileNotFoundError(f"Could not find {results_file}")

    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Load data
    df = pd.read_csv(results_file, sep="\t")
    df.columns = df.columns.str.strip()

    # Align pairwise change and pvalue columns 1:1 based on exact comparison prefixes
    pvalue_cols = [c for c in df.columns if c.endswith("_pvalue")]
    pairwise_prefixes = [c.replace("_pvalue", "") for c in pvalue_cols]
    change_cols = [f"{prefix}_change" for prefix in pairwise_prefixes]

    # Create a matching p-value matrix mapped to change column names
    pvalue_matrix = df[pvalue_cols].rename(
        columns=lambda c: c.replace("_pvalue", "_change")
    )

    # Mask change scores where the corresponding comparison p-value exceeds the threshold
    sig_change_matrix = df[change_cols].where(pvalue_matrix <= args.pvalue)

    # Keep only TFs with at least one statistically significant comparison
    valid_mask = sig_change_matrix.notna().any(axis=1)
    df = df[valid_mask].copy()
    sig_change_matrix = sig_change_matrix[valid_mask]

    if df.empty:
        raise ValueError(
            f"No TFs met the p-value threshold of <= {args.pvalue}"
        )

    # Rank TFs based strictly on maximum absolute change among significant comparisons
    df["max_abs_change"] = sig_change_matrix.abs().max(axis=1)

    # Select top N TFs with the largest significant differential changes
    top_tfs = df.sort_values(by="max_abs_change", ascending=False).head(args.n)

    # Set unique row index directly using output_prefix from BINDetect header
    heatmap_data = top_tfs.set_index("output_prefix")[change_cols]

    # Clean column headers for display
    heatmap_data.columns = [
        c.replace("_change", "") for c in heatmap_data.columns
    ]

    # Adjust figure height dynamically based on N
    fig_height = max(6, int(args.n * 0.4))
    plt.figure(figsize=(10, fig_height))

    # Plot actual differential binding change scores centered at 0
    sns.heatmap(
        heatmap_data,
        cmap="coolwarm",
        center=0,
        annot=True,
        fmt=".3f",
        linewidths=0.5,
        cbar_kws={"label": "Differential Binding Score (Change between samples)"},
    )
    plt.title(
        f"Top {args.n} Dynamic TFs: Differential Binding Change Between Samples (p <= {args.pvalue})",
        fontsize=12,
        fontweight="bold",
    )
    plt.ylabel("Transcription Factor (output_prefix)", fontsize=10)
    plt.tight_layout()

    out_heatmap_pdf = outdir / f"{args.prefix}_top{args.n}_differential_heatmap.pdf"
    out_heatmap_png = outdir / f"{args.prefix}_top{args.n}_differential_heatmap.png"

    plt.savefig(out_heatmap_pdf, dpi=300)
    plt.savefig(out_heatmap_png, dpi=300)
    plt.close()

    print(f"Saved: {out_heatmap_pdf}")
    print(f"Saved: {out_heatmap_png}")


if __name__ == "__main__":
    main()
