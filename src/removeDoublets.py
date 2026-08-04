#!/usr/bin/env python3

import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import snapatac2 as snap


def plot_hist(adata, column, prefix, stage, bins=50):
    if column not in adata.obs:
        print(f"Skipping {column}: not found")
        return
    plt.figure(figsize=(7, 5))
    vals = np.array(adata.obs[column])
    plt.hist(vals, bins=bins)
    plt.xlabel(column)
    plt.ylabel("Number of cells")
    plt.title(f"{stage}: {column}")
    plt.tight_layout()
    outfile = os.path.join("figures", f"{prefix}_{stage}_{column}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved {outfile}")


def plot_scatter(adata, x, y, prefix, stage):
    if x not in adata.obs or y not in adata.obs:
        return
    plt.figure(figsize=(7, 6))
    x_vals = np.array(adata.obs[x])
    y_vals = np.array(adata.obs[y])
    plt.scatter(x_vals, y_vals, s=1, alpha=0.3)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(f"{stage}: {y} vs {x}")
    plt.tight_layout()
    outfile = os.path.join("figures", f"{prefix}_{stage}_{y}_vs_{x}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description="Doublet removal using SnapATAC2's Scrublet implementation."
    )
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--output", required=True, help="Output h5ad file after doublet removal")
    parser.add_argument("--prefix", required=True, help="Prefix for output figures")
    args = parser.parse_args()

    # Create the figures directory
    os.makedirs("figures", exist_ok=True)

    # SIMPLE FIX 1: Ensure the destination directory for the output file actually exists
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print("Reading object...")
    adata = snap.read(args.input)
    print(f"Loaded: {adata.n_obs} cells")

    print("Selecting variable features...")
    snap.pp.select_features(adata)

    print("Running doublet detection (snap.pp.scrublet)...")
    snap.pp.scrublet(adata)

    print("Generating before-removal figures...")
    plot_hist(adata, "doublet_probability", args.prefix, "before_doublet_removal")
    plot_scatter(adata, "n_fragment", "doublet_probability", args.prefix, "before_doublet_removal")

    print("Filtering doublets...")
    snap.pp.filter_doublets(adata, probability_threshold=0.5)
    print(f"Remaining cells after doublet removal: {adata.n_obs}")

    print("Generating after-removal figures...")
    plot_hist(adata, "doublet_probability", args.prefix, "after_doublet_removal")
    plot_scatter(adata, "n_fragment", "doublet_probability", args.prefix, "after_doublet_removal")

    # SIMPLE FIX 2: Break the active read lock on the input file before writing
    print("Loading filtered data into memory for saving...")
    adata = adata.to_memory()

    print(f"Saving filtered object to {args.output}")
    adata.write(args.output)
    print("Done.")


if __name__ == "__main__":
    main()

