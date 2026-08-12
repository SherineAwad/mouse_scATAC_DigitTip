#!/usr/bin/env python3

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import snapatac2 as snap


def parse_peak_regions(df):
    regions = []
    for feature in df["feature name"]:
        chrom, coords = feature.split(":")
        start, end = coords.split("-")
        regions.append(f"{chrom}:{start}-{end}")
    return regions


def convert_to_pandas(result):
    if isinstance(result, pd.DataFrame):
        return result
    if hasattr(result, "to_pandas"):
        return result.to_pandas()
    raise TypeError(f"Unsupported motif enrichment result type: {type(result)}")


def plot_motif_enrichment(csv_file, prefix, label, pval_cutoff, top_n=20):
    df = pd.read_csv(csv_file)

    if df.empty:
        print(f"No motif results to plot for {label}")
        return

    required_columns = ["name", "log2(fold change)", "adjusted p-value"]
    for column in required_columns:
        if column not in df.columns:
            print(f"Cannot plot {label}: missing column '{column}'")
            return

    df = df.dropna(
        subset=["name", "log2(fold change)", "adjusted p-value"]
    )

    if df.empty:
        print(f"No valid motif results to plot for {label}")
        return

    # Keep statistically significant motifs
    significant_hits = df[
        df["adjusted p-value"] < pval_cutoff
    ].copy()

    if significant_hits.empty:
        print(
            f"No significant motif results to plot for {label}"
        )
        return

    # Rank using both statistical significance and effect size
    significant_hits["ranking_score"] = (
        -np.log10(
            np.maximum(
                significant_hits["adjusted p-value"],
                np.finfo(float).tiny
            )
        )
        * np.abs(
            significant_hits["log2(fold change)"]
        )
    )

    half_n = top_n // 2

    enriched_df = (
        significant_hits[
            significant_hits["log2(fold change)"] > 0
        ]
        .sort_values(
            "ranking_score",
            ascending=False
        )
        .head(half_n)
    )

    depleted_df = (
        significant_hits[
            significant_hits["log2(fold change)"] < 0
        ]
        .sort_values(
            "ranking_score",
            ascending=False
        )
        .head(half_n)
    )

    df = pd.concat(
        [enriched_df, depleted_df]
    ).sort_values(
        "log2(fold change)",
        ascending=True
    )

    if df.empty:
        print(
            f"No valid enriched or depleted motif results "
            f"passed significance to plot for {label}"
        )
        return

    plt.figure(figsize=(10, 8))

    colors = [
        "crimson" if x > 0 else "dodgerblue"
        for x in df["log2(fold change)"]
    ]

    plt.barh(
        df["name"],
        df["log2(fold change)"],
        color=colors
    )

    plt.axvline(
        0,
        linewidth=1.0,
        color="black",
        linestyle="-"
    )

    plt.xlabel("log2(fold change)")
    plt.ylabel("Motif")

    plt.title(
        f"{label.capitalize()} Motif Enrichment "
        f"(Red = Enriched | Blue = Depleted)"
    )

    plt.tight_layout()

    outfile = os.path.join(
        "figures",
        f"{prefix}_{label}_motif_enrichment.png"
    )

    plt.savefig(
        outfile,
        dpi=300
    )

    plt.close()

    print(
        f"Saved: {outfile}"
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Perform SnapATAC2 motif enrichment on "
            "significant differential peaks."
        )
    )

    parser.add_argument(
        "--peaks",
        required=True,
        help="Filtered annotated significant peak CSV"
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output files"
    )

    parser.add_argument(
        "--pval-cutoff",
        type=float,
        default=0.05,
        help=(
            "Adjusted p-value cutoff for significance "
            "(default: 0.05)"
        )
    )

    parser.add_argument(
        "--logfc-cutoff",
        type=float,
        default=0.25,
        help=(
            "Log2 fold change cutoff for significance "
            "(default: 0.25)"
        )
    )

    args = parser.parse_args()

    os.makedirs(
        "figures",
        exist_ok=True
    )

    print(
        f"Reading peaks: {args.peaks}"
    )

    df = pd.read_csv(
        args.peaks
    )

    if "feature name" not in df.columns:
        raise ValueError(
            "Required column 'feature name' not found "
            "in input CSV."
        )

    if "log2(fold_change)" not in df.columns:
        raise ValueError(
            "Required column 'log2(fold_change)' "
            "not found in input CSV."
        )

    if "adjusted p-value" not in df.columns:
        raise ValueError(
            "Required column 'adjusted p-value' "
            "not found in input CSV."
        )

    all_peaks = df.copy()

    all_regions = parse_peak_regions(
        all_peaks
    )

    print(
        f"Total peaks (background universe): "
        f"{len(all_regions)}"
    )

    gain_df = df[
        (df["log2(fold_change)"] > args.logfc_cutoff)
        &
        (df["adjusted p-value"] < args.pval_cutoff)
    ].copy()

    gain_regions = parse_peak_regions(
        gain_df
    )

    print(
        f"Significant gain peaks: "
        f"{len(gain_regions)}"
    )

    target_loss_cutoff = -abs(
        args.logfc_cutoff
    )

    loss_df = df[
        (df["log2(fold_change)"] < target_loss_cutoff)
        &
        (df["adjusted p-value"] < args.pval_cutoff)
    ].copy()

    loss_regions = parse_peak_regions(
        loss_df
    )

    print(
        f"Significant loss peaks: "
        f"{len(loss_regions)}"
    )

    regions = {
        "gain": gain_regions,
        "loss": loss_regions
    }

    print()
    print(
        "Loading mouse transcription-factor motifs..."
    )

    motifs = snap.datasets.cis_bp(
        unique=True
    )

    print(
        f"Loaded {len(motifs)} motifs."
    )

    print()
    print(
        "Running motif enrichment..."
    )

    enrichment = snap.tl.motif_enrichment(
        motifs=motifs,
        regions=regions,
        genome_fasta=snap.genome.mm39,
        background=all_regions,
        method="binomial",
    )

    print(
        "Motif enrichment completed."
    )

    for label in [
        "gain",
        "loss"
    ]:

        if label not in enrichment:
            print(
                f"No enrichment result for {label}"
            )
            continue

        print()
        print(
            f"Saving {label} results..."
        )

        result = enrichment[label]

        result = convert_to_pandas(
            result
        )

        result = result.replace(
            [np.inf, -np.inf],
            np.nan
        )

        result = result.dropna(
            subset=[
                "log2(fold change)",
                "adjusted p-value"
            ]
        )

        if result.empty:
            print(
                f"No valid motifs found for {label}. "
                f"Skipping."
            )
            continue

        csv_file = (
            f"{args.prefix}_"
            f"{label}_motif_enrichment.csv"
        )

        result.to_csv(
            csv_file,
            index=False
        )

        print(
            f"Saved: {csv_file}"
        )

        print(
            f"Generating {label} plot..."
        )

        plot_motif_enrichment(
            csv_file,
            args.prefix,
            label,
            args.pval_cutoff
        )

    print()
    print(
        "Done."
    )


if __name__ == "__main__":
    main()
