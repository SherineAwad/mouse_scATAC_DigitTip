#!/usr/bin/env python3

import argparse
import anndata
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Print summary statistics of an h5ad file."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input h5ad file"
    )
    args = parser.parse_args()

    adata = anndata.read_h5ad(args.input)

    print(f"File: {args.input}")
    print(f"Number of cells (n_obs): {adata.n_obs}")
    print(f"Number of features (n_vars): {adata.n_vars}")

    # If QC columns exist, print summary stats
    qc_cols = ["n_fragment", "frac_dup", "tsse"]

    existing = [
        col for col in qc_cols
        if col in adata.obs.columns
    ]

    if existing:
        print("\nQC metrics summary:")

        for col in existing:
            vals = adata.obs[col].values

            print(f"  {col}:")
            print(f"    min: {np.min(vals):.3f}")
            print(f"    max: {np.max(vals):.3f}")
            print(f"    mean: {np.mean(vals):.3f}")
            print(f"    median: {np.median(vals):.3f}")

    else:
        print("\nNo QC metrics found in obs.")

    # QC summary by sample
    if "sample" in adata.obs.columns:

        print("\nQC metrics by sample:")

        for sample in adata.obs["sample"].unique():

            sample_data = adata.obs[
                adata.obs["sample"] == sample
            ]

            print(f"\n  Sample: {sample}")
            print(f"    Cells: {len(sample_data)}")

            for col in existing:
                vals = sample_data[col].values

                print(
                    f"    {col} median: "
                    f"{np.median(vals):.3f}"
                )

    else:
        print(
            "\nNo 'sample' column found in adata.obs. "
            "Skipping sample-level QC."
        )

    # Check for layers
    if adata.layers:
        print(
            f"\nLayers present: "
            f"{list(adata.layers.keys())}"
        )
    else:
        print("\nNo layers present.")

    # Check for embeddings
    if adata.obsm:
        print(
            f"\nObsm keys: "
            f"{list(adata.obsm.keys())}"
        )
    else:
        print("\nNo obsm embeddings.")

    # Check for uns metadata
    if adata.uns:
        print(
            f"\nUns keys: "
            f"{list(adata.uns.keys())}"
        )


if __name__ == "__main__":
    main()
