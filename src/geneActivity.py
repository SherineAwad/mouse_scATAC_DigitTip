
#!/usr/bin/env python3

import argparse
import os

import pandas as pd
import scanpy as sc
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
    # 1. Read gene-activity matrix
    # ------------------------------------------------------------

    print("Reading gene matrix...")
    adata = sc.read_h5ad(args.input)

    print(
        f"Gene matrix: "
        f"{adata.n_obs} cells x "
        f"{adata.n_vars} genes"
    )

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("The input matrix contains no cells or genes.")

    # Require unique gene names for unambiguous marker results
    if not adata.var_names.is_unique:
        raise ValueError(
            "Gene names are not unique. Resolve duplicate gene identifiers "
            "before running marker analysis."
        )

    # Make cell identifiers unique without changing row order
    if not adata.obs_names.is_unique:

        print("Making duplicate cell identifiers unique...")

        if "sample" in adata.obs.columns:
            sample = adata.obs["sample"].astype(str)
            adata.obs_names = (
                sample + "_" + adata.obs_names.astype(str)
            )

        if not adata.obs_names.is_unique:
            adata.obs_names_make_unique()

    # Use existing log1p values; do not normalize them again
    if "log1p" in adata.layers:

        print("Using existing 'log1p' gene-activity layer...")
        adata.X = adata.layers["log1p"].copy()

    else:

        print(
            "No 'log1p' layer found. "
            "Normalizing and log-transforming gene activity..."
        )

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # ------------------------------------------------------------
    # 2. Check Leiden clusters
    # ------------------------------------------------------------

    if "leiden" not in adata.obs.columns:
        raise ValueError(
            "Input file does not contain adata.obs['leiden']."
        )

    if adata.obs["leiden"].isna().any():
        raise ValueError("The Leiden column contains missing values.")

    if adata.obs["leiden"].nunique() < 2:
        raise ValueError("At least two Leiden clusters are required.")

    print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

    # ------------------------------------------------------------
    # 3. Rank differential gene activity
    # ------------------------------------------------------------

    print("Calculating differential gene activity...")

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        use_raw=False,
        corr_method="benjamini-hochberg"
    )

    # ------------------------------------------------------------
    # 4. Extract top 30 markers per cluster
    # ------------------------------------------------------------

    n_top = 30
    n_plot = 10

    ranked = adata.uns["rank_genes_groups"]
    groups = ranked["names"].dtype.names

    results = []
    plot_genes_by_cluster = {}

    for cluster in groups:

        genes = ranked["names"][cluster]
        scores = ranked["scores"][cluster]
        pvals = ranked["pvals"][cluster]
        pvals_adj = ranked["pvals_adj"][cluster]

        for rank in range(min(n_top, len(genes))):

            results.append(
                {
                    "cluster": cluster,
                    "rank": rank + 1,
                    "gene": genes[rank],
                    "score": scores[rank],
                    "pval": pvals[rank],
                    "pval_adj": pvals_adj[rank],
                }
            )

        # Preserve the top 10 genes for each cluster.
        # Shared genes may appear in multiple cluster groups.
        plot_genes_by_cluster[cluster] = list(
            dict.fromkeys(genes[:n_plot].tolist())
        )

    markers = pd.DataFrame(results)

    marker_file = f"{args.prefix}_top_markers.csv"
    markers.to_csv(marker_file, index=False)

    print(f"Saved: {marker_file}")

    # ------------------------------------------------------------
    # 5. Create dotplot
    # ------------------------------------------------------------

    os.makedirs("figures", exist_ok=True)

    plot_file = f"figures/{args.prefix}_marker_dotplot.png"

    print("Creating marker dotplot...")

    dotplot = sc.pl.dotplot(
        adata,
        var_names=plot_genes_by_cluster,
        groupby="leiden",
        standard_scale="var",
        dendrogram=False,
        use_raw=False,
        show=False,
        return_fig=True
    )

    dotplot.savefig(
        plot_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close("all")

    print(f"Saved: {plot_file}")
    print("Done.")


if __name__ == "__main__":
    main()
