#!/usr/bin/env python3

import argparse
import os
import anndata

import matplotlib.pyplot as plt
import numpy as np
import snapatac2 as snap


def plot_hist(adata, column, prefix, stage):

    if column not in adata.obs.columns:
        print(f"Skipping {column}: not found")
        return

    plt.figure(figsize=(7, 5))

    plt.hist(
        adata.obs[column].values,
        bins=100
    )

    plt.xlabel(column)
    plt.ylabel("Number of cells")
    plt.title(f"{stage}: {column}")

    plt.tight_layout()

    outfile = os.path.join(
        "figures",
        f"{prefix}_{stage}_{column}.png"
    )

    plt.savefig(
        outfile,
        dpi=300
    )

    plt.close()

    print(f"Saved {outfile}")


def plot_scatter(adata, x, y, prefix, stage):

    if x not in adata.obs.columns:
        return

    if y not in adata.obs.columns:
        return

    plt.figure(figsize=(7, 6))

    plt.scatter(
        adata.obs[x],
        adata.obs[y],
        s=1,
        alpha=0.3
    )

    plt.xlabel(x)
    plt.ylabel(y)

    plt.title(
        f"{stage}: {y} vs {x}"
    )

    plt.tight_layout()

    outfile = os.path.join(
        "figures",
        f"{prefix}_{stage}_{y}_vs_{x}.png"
    )

    plt.savefig(
        outfile,
        dpi=300
    )

    plt.close()

    print(f"Saved {outfile}")


def make_qc_plots(adata, prefix, stage):

    qc_columns = [
        "n_fragment",
        "frac_dup",
        "tsse"
    ]

    for col in qc_columns:
        plot_hist(
            adata,
            col,
            prefix,
            stage
        )

    plot_scatter(
        adata,
        "n_fragment",
        "frac_dup",
        prefix,
        stage
    )

    plot_scatter(
        adata,
        "n_fragment",
        "tsse",
        prefix,
        stage
    )


def main():

    parser = argparse.ArgumentParser(
        description="snapATAC2 QC filtering workflow"
    )

    parser.add_argument(
        "--input",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    parser.add_argument(
        "--prefix",
        required=True
    )

    parser.add_argument(
        "--min-fragments",
        type=int,
        default=1000
    )

    parser.add_argument(
        "--max-fragments",
        type=int,
        default=100000
    )

    parser.add_argument(
        "--min-tsse",
        type=float,
        default=1.5,
        help="Minimum TSS enrichment score (default: 1.5)"
    )

    parser.add_argument(
        "--max-tsse",
        type=float,
        default=100,
        help="Maximum TSS enrichment score (default: 100, effectively no upper bound)"
    )

    parser.add_argument(
        "--max-frac-dup",
        type=float,
        default=0.5,
        help="Maximum fraction of duplicate fragments (default: 0.5)"
    )

    args = parser.parse_args()

    os.makedirs(
        "figures",
        exist_ok=True
    )

    print("Reading object")

    adata = anndata.read_h5ad(args.input)

    print(adata)

    print(
        "Generating pre-filter QC plots"
    )

    make_qc_plots(
        adata,
        args.prefix,
        "prefilter"
    )

    print(
        "Filtering cells"
    )

    keep = (
        (adata.obs["n_fragment"] >= args.min_fragments)
        &
        (adata.obs["n_fragment"] <= args.max_fragments)
        &
        (adata.obs["tsse"] >= args.min_tsse)
        &
        (adata.obs["tsse"] <= args.max_tsse)
        &
        (adata.obs["frac_dup"] <= args.max_frac_dup)
    )

    print(
        f"Keeping {keep.sum()} / {adata.n_obs} cells"
    )

    print(
        "Subsetting cells"
    )

    adata = adata[keep].copy()

    print(
        "Generating post-filter QC plots"
    )

    make_qc_plots(
        adata,
        args.prefix,
        "postfilter"
    )

    print(
        "Adding tile matrix after filtering"
    )

    snap.pp.add_tile_matrix(
        adata
    )

    print(
        "Writing filtered object"
    )

    adata.write(
        args.output,
        compression="gzip"
    )

    print("Done")


if __name__ == "__main__":
    main()
