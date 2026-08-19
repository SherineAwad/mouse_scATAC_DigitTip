import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import snapatac2 as snap


def main():

    parser = argparse.ArgumentParser(
        description="Differential accessibility analysis on called peaks (MACS3)."
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--group1", required=True)
    parser.add_argument("--group2", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--pval-cutoff", type=float, default=0.05)
    parser.add_argument("--logfc-cutoff", type=float, default=0.5)

    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    print("Reading object...")
    adata = snap.read(args.input)

    print(f"Loaded: {adata.n_obs} cells")

    # ------------------------------------------------------------
    # 1. Call peaks separately for each sample
    # ------------------------------------------------------------

    print("Calling MACS3 peaks per sample...")

    peaks = snap.tl.macs3(
        adata,
        groupby="sample",
        qvalue=0.05,
        n_jobs=8,
        inplace=False
    )

    print("MACS3 peak calling complete.")

    for sample, peak_df in peaks.items():
        print(
            f"  {sample}: {peak_df.height} peaks"
        )

    # ------------------------------------------------------------
    # 2. Merge peaks into a common peak set
    # ------------------------------------------------------------

    print("Merging peaks...")

    merged_peaks = snap.tl.merge_peaks(
        peaks,
        snap.genome.mm39
    )

    print(
        f"Common peak set: {merged_peaks.height} peaks"
    )

    # ------------------------------------------------------------
    # 3. Save merged peaks as BED
    # ------------------------------------------------------------

    bed_file = f"{args.prefix}_merged_peaks.bed"

    with open(bed_file, "w") as f:

        for peak in merged_peaks["Peaks"].to_list():

            chrom, coordinates = peak.split(":")
            start, end = coordinates.split("-")

            f.write(
                f"{chrom}\t{start}\t{end}\n"
            )

    print(f"Saved {bed_file}")

    # ------------------------------------------------------------
    # 4. Create peak matrix
    # ------------------------------------------------------------

    print("Building peak matrix...")

    peak_mat = snap.pp.make_peak_matrix(
        adata,
        use_rep=merged_peaks["Peaks"].to_list(),
        inplace=False
    )

    print(
        f"Peak matrix: "
        f"{peak_mat.n_obs} cells x "
        f"{peak_mat.n_vars} peaks"
    )

    # ------------------------------------------------------------
    # 5. Get group indices
    # ------------------------------------------------------------

    group1_idx = np.where(
        adata.obs["sample"].to_numpy() == args.group1
    )[0]

    group2_idx = np.where(
        adata.obs["sample"].to_numpy() == args.group2
    )[0]

    if len(group1_idx) == 0:
        raise ValueError(
            f"Group1 '{args.group1}' has no cells."
        )

    if len(group2_idx) == 0:
        raise ValueError(
            f"Group2 '{args.group2}' has no cells."
        )

    print(
        f"Group1 ({args.group1}): {len(group1_idx)} cells"
    )

    print(
        f"Group2 ({args.group2}): {len(group2_idx)} cells"
    )

    # ------------------------------------------------------------
    # 6. Differential accessibility test
    # ------------------------------------------------------------

    print("Running differential accessibility test...")

    res = snap.tl.diff_test(
        peak_mat,
        cell_group1=group1_idx,
        cell_group2=group2_idx,
        features=None,
        direction="both",
        min_log_fc=0.0,
        min_pct=0.05
    )

    print(
        f"Tested {res.height} peaks."
    )

    # ------------------------------------------------------------
    # 7. Save all results
    # ------------------------------------------------------------

    all_results = (
        f"{args.prefix}_all_results.csv"
    )

    res.write_csv(all_results)

    print(
        f"Saved {all_results}"
    )

    # ------------------------------------------------------------
    # 8. Significant peaks
    # ------------------------------------------------------------

    pval_col = "adjusted p-value"
    logfc_col = "log2(fold_change)"

    res_sig = res.filter(
        (res[pval_col] < args.pval_cutoff)
        &
        (
            np.abs(res[logfc_col])
            >= args.logfc_cutoff
        )
    )

    up = res_sig.filter(
        res_sig[logfc_col] > 0
    )

    down = res_sig.filter(
        res_sig[logfc_col] < 0
    )

    significant_results = (
        f"{args.prefix}_significant_results.csv"
    )

    res_sig.write_csv(
        significant_results
    )

    print(
        f"Saved {significant_results}"
    )

    print(
        f"Significant peaks: {res_sig.height}"
    )

    print(
        f"Higher in {args.group1}: {up.height}"
    )

    print(
        f"Higher in {args.group2}: {down.height}"
    )

    # ------------------------------------------------------------
    # 9. Volcano plot
    # ------------------------------------------------------------

    print("Generating volcano plot...")

    pvalues = np.maximum(
        res[pval_col].to_numpy(),
        np.finfo(float).tiny
    )

    plt.figure(figsize=(8, 6))

    plt.scatter(
        res[logfc_col],
        -np.log10(pvalues),
        s=2,
        alpha=0.3,
        color="gray"
    )

    if up.height > 0:
        up_pvalues = np.maximum(
            up[pval_col].to_numpy(),
            np.finfo(float).tiny
        )

        plt.scatter(
            up[logfc_col],
            -np.log10(up_pvalues),
            s=2,
            alpha=0.7,
            color="red",
            label=f"Higher in {args.group1}"
        )

    if down.height > 0:
        down_pvalues = np.maximum(
            down[pval_col].to_numpy(),
            np.finfo(float).tiny
        )

        plt.scatter(
            down[logfc_col],
            -np.log10(down_pvalues),
            s=2,
            alpha=0.7,
            color="blue",
            label=f"Higher in {args.group2}"
        )

    plt.axhline(
        -np.log10(args.pval_cutoff),
        linestyle="--",
        color="black",
        alpha=0.5
    )

    plt.axvline(
        args.logfc_cutoff,
        linestyle="--",
        color="black",
        alpha=0.5
    )

    plt.axvline(
        -args.logfc_cutoff,
        linestyle="--",
        color="black",
        alpha=0.5
    )

    plt.xlabel("log2 fold change")
    plt.ylabel("-log10 adjusted p-value")

    plt.title(
        f"Differential accessibility: "
        f"{args.group1} vs {args.group2}"
    )

    plt.legend()

    plt.tight_layout()

    outfile = os.path.join(
        "figures",
        f"{args.prefix}_volcano.png"
    )

    plt.savefig(
        outfile,
        dpi=300
    )

    plt.close()

    print(
        f"Saved {outfile}"
    )

    # ------------------------------------------------------------
    # 10. Save peak matrix
    # ------------------------------------------------------------

    print(
        f"Saving peak matrix to {args.output}"
    )

    contrast_cells = np.concatenate([group1_idx, group2_idx])
    peak_mat_out = peak_mat[contrast_cells, :].copy()

    peak_mat_out.write(args.output)

    print("Done.")


if __name__ == "__main__":
    main()
