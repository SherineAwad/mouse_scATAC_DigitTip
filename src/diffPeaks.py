#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import snapatac2 as snap

def main():
    parser = argparse.ArgumentParser(
        description="Differential accessibility analysis between two cell groups."
    )
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--output", required=True, help="Output h5ad file with diff test results")
    parser.add_argument("--group1", required=True, help="Group label (sample name or cluster)")
    parser.add_argument("--group2", required=True, help="Group label (sample name or cluster)")
    parser.add_argument("--prefix", required=True, help="Prefix for output figures")
    parser.add_argument("--pval-cutoff", type=float, default=0.05, help="Adjusted p-value cutoff (default: 0.05)")
    parser.add_argument("--logfc-cutoff", type=float, default=0.5, help="Log2 fold change cutoff (default: 0.5)")
    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    print("Reading object...")
    adata = snap.read(args.input)
    print(f"Loaded: {adata.n_obs} cells")

    # Get indices for each group (using 'sample' column)
    group1_idx = np.where(adata.obs["sample"] == args.group1)[0]
    group2_idx = np.where(adata.obs["sample"] == args.group2)[0]

    if len(group1_idx) == 0:
        raise ValueError(f"Group1 '{args.group1}' has no cells.")
    if len(group2_idx) == 0:
        raise ValueError(f"Group2 '{args.group2}' has no cells.")

    print(f"Group1 ({args.group1}): {len(group1_idx)} cells")
    print(f"Group2 ({args.group2}): {len(group2_idx)} cells")

    # Run differential test
    print("Running differential accessibility test...")
    res = snap.tl.diff_test(
        adata,
        cell_group1=group1_idx,
        cell_group2=group2_idx,
        features=None,
        direction="both",
        min_log_fc=args.logfc_cutoff,
        min_pct=0.05
        # solver removed
    )

    # Store results in uns
    adata.uns["diffPeaks"] = res

    # Generate volcano plot
    print("Generating volcano plot...")
    pval_cut = args.pval_cutoff
    logfc_cut = args.logfc_cutoff

    # Use correct column names: 'pvalue_adj' and 'log2_fc'
    res_sig = res[(res["pvalue_adj"] < pval_cut) & (np.abs(res["log2_fc"]) > logfc_cut)]
    up = res_sig[res_sig["log2_fc"] > 0]
    down = res_sig[res_sig["log2_fc"] < 0]

    plt.figure(figsize=(8, 6))
    plt.scatter(res["log2_fc"], -np.log10(res["pvalue_adj"]), s=2, alpha=0.3, color="gray")
    plt.scatter(up["log2_fc"], -np.log10(up["pvalue_adj"]), s=2, alpha=0.7, color="red", label=f"Group1 ({args.group1})")
    plt.scatter(down["log2_fc"], -np.log10(down["pvalue_adj"]), s=2, alpha=0.7, color="blue", label=f"Group2 ({args.group2})")
    plt.axhline(y=-np.log10(pval_cut), linestyle="--", color="black", alpha=0.5, label=f"p‑adj = {pval_cut}")
    plt.axvline(x=logfc_cut, linestyle="--", color="black", alpha=0.5)
    plt.axvline(x=-logfc_cut, linestyle="--", color="black", alpha=0.5)
    plt.xlabel("log2 fold change")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(f"Differential accessibility: {args.group1} vs {args.group2}")
    plt.legend()
    plt.tight_layout()
    outfile = os.path.join("figures", f"{args.prefix}_volcano.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved {outfile}")

    # Save final object
    print(f"Saving object to {args.output}")
    adata.write(args.output)
    print("Done.")

if __name__ == "__main__":
    main()
