#!/usr/bin/env python3

import argparse
import os
import anndata
import snapatac2 as snap
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description="Spectral embedding, clustering, UMAP, and tile matrix layers."
    )
    parser.add_argument("--input", required=True, help="Input h5ad file (after doublet removal)")
    parser.add_argument("--output", required=True, help="Output h5ad file with embeddings, clusters, and tile layers")
    parser.add_argument("--prefix", required=True, help="Prefix for output figures")
    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    print("Reading object...")
    adata = anndata.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells")

    # Ensure we have a tile matrix (adata.X should be present after add_tile_matrix)
    if adata.X is None:
        raise RuntimeError("No count matrix found. Did you run snap.pp.add_tile_matrix?")

    # ---- Store raw counts as layer ----
    print("Storing raw tile counts as layer 'counts'...")
    adata.layers["counts"] = adata.X.copy()

    # ---- Log1p normalization on the tile matrix ----
    print("Computing log1p and storing as layer 'log1p'...")
    # Copy to avoid modifying original X
    X_log1p = adata.X.copy()
    # Normalize by total per cell (like sc.pp.normalize_total) and log1p
    # We'll use scanpy functions on a temporary AnnData
    tmp = sc.AnnData(X=X_log1p)
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    adata.layers["log1p"] = tmp.X.copy()
    del tmp

    # ---- Now proceed with clustering on the original tile matrix (X) ----
    print("Running spectral embedding...")
    snap.tl.spectral(adata)

    print("Computing nearest-neighbor graph...")
    snap.pp.knn(adata, use_rep="X_spectral")

    print("Running Leiden clustering...")
    snap.tl.leiden(adata)

    print("Running UMAP...")
    snap.tl.umap(adata)

    print("Generating UMAP figure...")
    plt.figure(figsize=(18, 8))
    sc.pl.umap(adata, color="leiden", show=False)
    outfile = os.path.join("figures", f"{args.prefix}_umap_clusters.png")
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {outfile}")
    
    print("Generating QC UMAP plots...")

    sc.pl.umap(
        adata,
        color=[
        "leiden",
        "sample",
        "n_fragment",
        "tsse",
        "frac_dup"
        ],
       ncols=3,
       show=False)

    outfile = os.path.join("figures", f"{args.prefix}_umap_qc.png")
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved {outfile}")

    print(f"Saving object to {args.output}")
    adata.write(args.output, compression="gzip")
    print("Done.")

    # QC metrics by Leiden cluster
    # -----------------------------

    if "leiden" in adata.obs.columns:

        qc_metrics = [
           "n_fragment",
           "tsse",
           "frac_dup"]

        fig, axes = plt.subplots(
            1,
            len(qc_metrics),
            figsize=(18, 5)
        )

        for ax, metric in zip(axes, qc_metrics):

            if metric in adata.obs.columns:

               (
                    adata.obs
                    .groupby("leiden")[metric]
                    .median()
                    .sort_index()
                    .plot(
                       kind="bar",
                       ax=ax
                    )
                )

               ax.set_title(f"Median {metric} per Leiden cluster")
               ax.set_xlabel("Leiden cluster")
               ax.set_ylabel(metric)
               ax.tick_params(axis="x", rotation=90)

        plt.tight_layout()

        outfile = os.path.join(
          "figures",
          f"{args.prefix}_qc_metrics_by_leiden.png"
         )

        plt.savefig(
         outfile,
         dpi=300,
         bbox_inches="tight"
        )

        plt.close()

        print(f"Saved {outfile}")

if __name__ == "__main__":
    main()
