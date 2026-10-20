import argparse
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    parser = argparse.ArgumentParser(
        description="Plot absolute mean TF activity per sample (capturing activated and repressed TFs) with dotplot and clustermap."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input chromVAR.h5ad file"
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output figures and summary CSV"
    )
    parser.add_argument(
        "--n-top",
        type=int,
        default=5,
        help="Number of top absolute TFs to select per sample (default: 5)"
    )
    
    args = parser.parse_args()

    # Ensure figures directory exists
    os.makedirs("figures", exist_ok=True)

    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    adata.obs_names_make_unique()

    # Validate 'sample' column exists
    if "sample" not in adata.obs.columns:
        raise ValueError("The 'sample' column is missing from adata.obs.")

    # Enforce strict sample order for Scanpy and Pandas
    desired_order = ["control", "1xAmp", "3xAmp", "5xAmp"]
    adata.obs['sample'] = pd.Categorical(
        adata.obs['sample'], 
        categories=desired_order, 
        ordered=True
    )

    print("Calculating mean motif activity per sample...")
    
    # Extract matrix and group by sample to compute mean activity
    X_mat = adata.X.toarray() if hasattr(adata.X, "toarray") else np.array(adata.X)
    expr_df = pd.DataFrame(X_mat, index=adata.obs_names, columns=adata.var_names)
    expr_df['sample'] = adata.obs['sample'].values

    mean_df = expr_df.groupby('sample', observed=False).mean()
    mean_df = mean_df.reindex(desired_order)

    # Find top N TFs independently based on absolute magnitude (activated & repressed) and take their union
    selected_tfs = set()
    for sample_name, row in mean_df.iterrows():
        top_sample_tfs = row.abs().nlargest(args.n_top).index.tolist()
        selected_tfs.update(top_sample_tfs)

    selected_tfs = list(selected_tfs)
    print(f"Identified {len(selected_tfs)} unique top absolute TFs across all samples.")

    # Save the mean activity table for these selected TFs
    csv_output = f"{args.prefix}_mean_activity_per_sample.csv"
    mean_df.T.loc[selected_tfs].to_csv(csv_output)
    print(f"Saved mean activity summary table to {csv_output}")

    # 1. Generate and save the Scanpy dotplot
    dotplot_output = os.path.join("figures", f"{args.prefix}_dotplot.png")
    sc.pl.dotplot(
        adata,
        var_names=selected_tfs,
        groupby="sample",
        dendrogram=False,
        show=False
    )
    plt.savefig(dotplot_output, bbox_inches="tight")
    plt.close()
    print(f"Saved dotplot figure to {dotplot_output}")

    # 2. Generate and save the Seaborn clustermap for absolute/divergent Z-scores
    clustermap_output = os.path.join("figures", f"{args.prefix}_clustermap.png")
    data_to_plot = mean_df.T.loc[selected_tfs]
    
    sns.clustermap(
        data_to_plot,
        cmap="coolwarm",
        center=0,
        col_cluster=False,
        figsize=(10, max(6, len(selected_tfs) * 0.25)),
        dendrogram_ratio=0.15,
        cbar_kws={"label": "Mean ChromVAR Deviation"}
    )
    plt.savefig(clustermap_output, bbox_inches="tight")
    plt.close()
    print(f"Saved clustermap figure to {clustermap_output}")
    
    print("Done!")

if __name__ == "__main__":
    main()
