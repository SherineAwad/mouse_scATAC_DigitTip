#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import matplotlib.pyplot as plt


def main():

    parser = argparse.ArgumentParser(
        description="Find top differential gene-activity markers by Leiden cluster."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input gene activity H5AD file"
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output files"
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # 1. Read gene matrix
    # ------------------------------------------------------------

    print("Reading gene matrix...")
    adata = snap.read(args.input)
    adata = adata.to_memory()

    print(
        f"Gene matrix: "
        f"{adata.n_obs} cells x "
        f"{adata.n_vars} genes"
    )

    # Use existing log-normalized gene activity values if available.
    # Otherwise create them from the gene activity matrix in .X.
    if "log1p" in adata.layers:

        print("Using existing 'log1p' gene-activity layer...")
        adata.X = adata.layers["log1p"].copy()

    else:

        print(
            "No 'log1p' layer found. "
            "Normalizing gene activity and computing log1p..."
        )

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    # ------------------------------------------------------------
    # 2. Check Leiden clusters
    # ------------------------------------------------------------

    if "leiden" not in adata.obs.columns:
        raise ValueError(
            "Input file does not contain adata.obs['leiden']"
        )

    print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

    # ------------------------------------------------------------
    # 3. Differential gene activity
    # ------------------------------------------------------------

    print("Calculating differential gene activity...")

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon"
    )

    # ------------------------------------------------------------
    # 4. Extract top markers
    # ------------------------------------------------------------

    n_top = 30

    results = []

    groups = adata.uns["rank_genes_groups"]["names"].dtype.names

    for cluster in groups:

        genes = adata.uns["rank_genes_groups"]["names"][cluster][:n_top]
        scores = adata.uns["rank_genes_groups"]["scores"][cluster][:n_top]
        logfc = adata.uns["rank_genes_groups"]["logfoldchanges"][cluster][:n_top]
        pvals = adata.uns["rank_genes_groups"]["pvals"][cluster][:n_top]
        pvals_adj = adata.uns["rank_genes_groups"]["pvals_adj"][cluster][:n_top]

        for rank, (gene, score, lf, pval, pval_adj) in enumerate(
            zip(genes, scores, logfc, pvals, pvals_adj),
            start=1
        ):
            results.append(
                {
                    "cluster": cluster,
                    "rank": rank,
                    "gene": gene,
                    "score": score,
                    "logfoldchange": lf,
                    "pval": pval,
                    "pval_adj": pval_adj,
                }
            )

    markers = pd.DataFrame(results)

    marker_file = f"{args.prefix}_top_markers.csv"
    markers.to_csv(marker_file, index=False)

    print(f"Saved: {marker_file}")

    # ------------------------------------------------------------
    # 5. Select genes for dot plot
    # ------------------------------------------------------------

    plot_genes = (
        markers
        .groupby("cluster", sort=False)
        .head(10)["gene"]
        .drop_duplicates()
        .tolist()
    )

    print(f"Genes in dot plot: {len(plot_genes)}")

    # ------------------------------------------------------------
    # 6. Dot plot
    # ------------------------------------------------------------

    os.makedirs("figures", exist_ok=True)

    print("Creating marker dot plot...")

    sc.pl.dotplot(
        adata,
        var_names=plot_genes,
        groupby="leiden",
        standard_scale="var",
        dendrogram=False,
        show=False
    )

    plt.savefig(
        f"figures/{args.prefix}_marker_dotplot.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(
        f"Saved: figures/{args.prefix}_marker_dotplot.png"
    )

    print("Done.")


if __name__ == "__main__":
    main()
