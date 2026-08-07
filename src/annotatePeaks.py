#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import snapatac2 as snap
import pyranges as pr

def main():
    parser = argparse.ArgumentParser(
        description="Annotate peaks with nearest genes using a GTF file."
    )
    parser.add_argument("--peaks", required=True, help="CSV file with peaks (must contain 'feature name' column)")
    parser.add_argument("--prefix", required=True, help="Prefix for output files")
    parser.add_argument("--dist-cutoff", type=int, default=50000, help="Distance cutoff (bp) for gene assignment (default: 50000)")
    parser.add_argument("--n", type=int, default=20, help="Number of top genes to label on volcano plot (default: 20)")
    args = parser.parse_args()

    # Load peaks from CSV
    print(f"Reading peaks from {args.peaks}...")
    peaks_df = pd.read_csv(args.peaks)
    if "feature name" not in peaks_df.columns:
        raise ValueError("CSV must contain a column named 'feature name' with peak coordinates.")
    
    peak_coords = peaks_df["feature name"].str.strip().tolist()

    # Get GFF3 from snapATAC2 mm39
    gff3_path = snap.genome.mm39.annotation

    if not os.path.exists(gff3_path):
       raise FileNotFoundError(f"GFF3 file not found: {gff3_path}")

    print(f"Reading GFF3: {gff3_path}...")
    gtf = pr.read_gff3(gff3_path)

    
    genes = gtf[gtf.Feature == "gene"].df.copy()
    # Ensure required columns exist
    print(genes.columns)
    genes["gene_name"] = genes["gene_name"] if "gene_name" in genes.columns else genes["gene_id"]

    # Convert peaks to PyRanges object
    peaks_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": [x.split(":")[0] for x in peak_coords],
        "Start": [int(x.split(":")[1].split("-")[0]) for x in peak_coords],
        "End": [int(x.split(":")[1].split("-")[1]) for x in peak_coords],
        "Peak": peak_coords
    }))

    # Find nearest gene within distance cutoff (TSS-based)
    print(f"Annotating peaks (distance cutoff: {args.dist_cutoff} bp)...")
    # For each peak, find the nearest TSS (using gene start for + strand, end for - strand)
    # We'll use nearest TSS distance via PyRanges' nearest() method on TSS positions.
    # Find nearest gene within distance cutoff (TSS-based)
    tss = genes.copy()
    tss["Start"] = np.where(tss["Strand"] == "+", tss["Start"], tss["End"])
    tss["End"] = tss["Start"]  # TSS is a single point

    # Compute distance to nearest TSS
    # Use PyRanges' nearest() with distance calculation
    tss = pr.PyRanges(tss)
    annotated = peaks_pr.nearest(tss, how="inner", overlap=False)
    # annotated has columns: Peak, gene_name, distance (if computed)
    # distance is from peak center to TSS? We'll compute actual distance manually later.

    # Now we need to manually compute distance from peak center to TSS
    # Since PyRanges nearest() returns distance between intervals (start/end), we'll use it.
    # However, the distance column is not automatically added; we'll compute it.
    # A simpler approach: use pyranges' distance to TSS with a custom function.

    # We'll use a loop (inefficient but fine for up to 100k peaks).
    print("Computing distances to nearest gene TSS...")
    results = []
    for _, peak_row in peaks_pr.df.iterrows():
        peak_chr = peak_row["Chromosome"]
        peak_center = (peak_row["Start"] + peak_row["End"]) // 2
        # Subset genes on same chromosome
        genes_chr = genes[genes["Chromosome"] == peak_chr]
        if genes_chr.empty:
            results.append({"Peak": peak_row["Peak"], "gene_name": None, "distance": None})
            continue
        # Compute TSS position
        tss_pos = np.where(genes_chr["Strand"] == "+", genes_chr["Start"], genes_chr["End"])
        dists = np.abs(tss_pos - peak_center)
        min_idx = np.argmin(dists)
        if dists[min_idx] <= args.dist_cutoff:
            results.append({
                "Peak": peak_row["Peak"],
                "gene_name": genes_chr["gene_name"].values[min_idx],
                "distance": int(dists[min_idx])
            })
        else:
            results.append({"Peak": peak_row["Peak"], "gene_name": None, "distance": None})

    results_df = pd.DataFrame(results)

    # Merge with original peaks_df to keep all columns
    final_df = peaks_df.merge(results_df, left_on="feature name", right_on="Peak", how="left")
    final_df.drop(columns=["Peak"], inplace=True, errors="ignore")

    # Save output
    outfile = f"{args.prefix}_annotated_peaks.csv"
    print(f"Saving annotated peaks to {outfile}...")
    final_df.to_csv(outfile, index=False)

    # ---- Volcano plot with gene labels (added) ----
    print("Generating volcano plot with gene labels...")
    logfc_col = "log2(fold_change)"
    pval_col = "adjusted p-value"

    if logfc_col in final_df.columns and pval_col in final_df.columns:
        # ---- ONLY CHANGE: filter to keep rows with gene_name ----
        final_df_with_gene = final_df[final_df["gene_name"].notna()].copy()
        final_df_sorted = final_df_with_gene.sort_values(pval_col, ascending=True).copy()
        top_n = final_df_sorted.head(args.n)

        plt.figure(figsize=(10, 8))
        plt.scatter(final_df[logfc_col], -np.log10(final_df[pval_col]), s=5, alpha=0.4, color="gray")
        plt.scatter(top_n[logfc_col], -np.log10(top_n[pval_col]), s=20, alpha=0.8, color="red")

        for _, row in top_n.iterrows():
            # gene_name is guaranteed to exist here
            label = row["gene_name"]
            plt.text(row[logfc_col], -np.log10(row[pval_col]), label, fontsize=6, alpha=0.8)

        plt.axhline(y=-np.log10(0.05), linestyle="--", color="black", alpha=0.5, label="p‑adj = 0.05")
        plt.axvline(x=0.5, linestyle="--", color="black", alpha=0.5)
        plt.axvline(x=-0.5, linestyle="--", color="black", alpha=0.5)
        plt.xlabel("log2 fold change")
        plt.ylabel("-log10 adjusted p-value")
        plot_title = f"{args.prefix}: Volcano plot with top {args.n} gene labels"
        plt.title(plot_title)
        plt.legend()
        plt.tight_layout()

        volc_out = os.path.join("figures", f"{args.prefix}_volcano_annotated.png")
        plt.savefig(volc_out, dpi=300)
        plt.close()
        print(f"Saved volcano plot to {volc_out}")
    else:
        print("Skipping volcano plot: required columns not found.")

    print("Done.")
    print(final_df["gene_name"].isna().sum())
    print(final_df["gene_name"].notna().sum()) 
    print("Distance cutoff:", args.dist_cutoff)

if __name__ == "__main__":
    main()
