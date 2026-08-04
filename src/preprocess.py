#!/usr/bin/env python3

import argparse
import gc
import gzip
import os
from collections import Counter, defaultdict
import bisect

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import snapatac2 as snap


def read_samples(sample_file):
    with open(sample_file) as f:
        return [line.strip() for line in f if line.strip()]


# ================== ORIGINAL PLOTTING FUNCTIONS (KEPT BUT NOT CALLED) ==================
def plot_hist(adata, column, prefix, stage):
    if column not in adata.obs.columns:
        raise RuntimeError(f"Column '{column}' missing for {prefix}")
    plt.figure(figsize=(7, 5))
    plt.hist(adata.obs[column].values, bins=100)
    plt.xlabel(column)
    plt.ylabel("Number of cells")
    plt.title(f"{prefix} {stage}: {column}")
    plt.tight_layout()
    outfile = os.path.join("figures", f"{prefix}_{stage}_{column}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved {outfile}")


def plot_scatter(adata, x, y, prefix, stage):
    if x not in adata.obs.columns or y not in adata.obs.columns:
        raise RuntimeError(f"Columns '{x}' or '{y}' missing for {prefix}")
    plt.figure(figsize=(7, 6))
    plt.scatter(adata.obs[x], adata.obs[y], s=1, alpha=0.3)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(f"{prefix} {stage}: {y} vs {x}")
    plt.tight_layout()
    outfile = os.path.join("figures", f"{prefix}_{stage}_{y}_vs_{x}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved {outfile}")


# ================== FRAGMENT SIZE DISTRIBUTION (UNIQUE) ==================
def compute_fragment_size_distribution(fragment_file):
    size_counter = Counter()
    with gzip.open(fragment_file, 'rt') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            start = int(parts[1])
            end = int(parts[2])
            length = end - start
            size_counter[length] += 1
    return size_counter


def plot_frag_size_dist_from_counter(size_counter, prefix, stage):
    sizes = sorted(size_counter.keys())
    counts = [size_counter[s] for s in sizes]
    plt.figure(figsize=(10, 5))
    plt.plot(sizes, counts, linewidth=1, color='black')
    plt.axvspan(80, 300, alpha=0.2, color='green', label='NFR (80–300 bp)')
    plt.axvspan(147, 200, alpha=0.2, color='blue', label='Mono‑nucleosome (147–200 bp)')
    plt.axvspan(300, 400, alpha=0.2, color='red', label='Di‑nucleosome (300–400 bp)')
    plt.xlabel("Fragment size (bp)")
    plt.ylabel("Number of fragments")
    plt.title(f"{prefix} {stage}: Fragment size distribution")
    plt.legend()
    plt.tight_layout()
    outfile = os.path.join("figures", f"{prefix}_{stage}_fragment_size.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved {outfile}")


# ================== TSS PROFILE (UNIQUE) ==================
def compute_tss_profile(fragment_file, genome, window=2000, bin_size=10):
    gtf_path = genome.annotation
    genes = pd.read_csv(gtf_path, sep='\t', comment='#', header=None,
                        names=['seqname', 'source', 'feature', 'start', 'end',
                               'score', 'strand', 'frame', 'attributes'])
    genes = genes[genes['feature'] == 'gene']

    tss_positions = []
    for idx, row in genes.iterrows():
        chrom = row['seqname']
        if row['strand'] == '+':
            tss = row['start']
        else:
            tss = row['end']
        tss_positions.append((chrom, tss))

    tss_by_chrom = defaultdict(list)
    for chrom, pos in tss_positions:
        tss_by_chrom[chrom].append(pos)
    for chrom in tss_by_chrom:
        tss_by_chrom[chrom].sort()

    n_bins = (2 * window) // bin_size + 1
    profile = np.zeros(n_bins)
    total_tss = 0

    with gzip.open(fragment_file, 'rt') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom = parts[0]
            if chrom not in tss_by_chrom:
                continue
            start = int(parts[1])
            end = int(parts[2])
            mid = (start + end) // 2

            tss_list = tss_by_chrom[chrom]
            i = bisect.bisect_left(tss_list, mid)
            nearest_dist = float('inf')
            if i < len(tss_list):
                dist = abs(tss_list[i] - mid)
                if dist < nearest_dist:
                    nearest_dist = dist
            if i > 0:
                dist = abs(tss_list[i-1] - mid)
                if dist < nearest_dist:
                    nearest_dist = dist

            if nearest_dist <= window:
                bin_idx = int((nearest_dist + window) // bin_size)
                if bin_idx < n_bins:
                    profile[bin_idx] += 1
                    total_tss += 1

    if total_tss > 0:
        profile = profile / profile.max()

    positions = np.arange(-window, window + bin_size, bin_size)[:n_bins]
    return positions, profile


def plot_tss_profile(positions, profile, prefix, stage):
    plt.figure(figsize=(10, 5))
    plt.plot(positions, profile, linewidth=1, color='black')
    plt.axvline(0, color='red', linestyle='--', linewidth=1, label='TSS')
    plt.xlabel("Distance to TSS (bp)")
    plt.ylabel("Relative coverage")
    plt.title(f"{prefix} {stage}: TSS enrichment profile")
    plt.legend()
    plt.tight_layout()
    outfile = os.path.join("figures", f"{prefix}_{stage}_tss_profile.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description="Import Cell Ranger scATAC fragments into one h5ad file per sample."
    )
    parser.add_argument(
        "--samples",
        required=True,
        help="Text file containing one sample name per line."
    )
    args = parser.parse_args()

    samples = read_samples(args.samples)
    os.makedirs("figures", exist_ok=True)

    for sample in samples:
        fragments = os.path.join(sample, "outs", "fragments.tsv.gz")
        if not os.path.exists(fragments):
            raise FileNotFoundError(f"Cannot find {fragments}")

        outfile = f"{sample}.h5ad"

        print("=" * 60)
        print(f"Processing {sample}")
        print(f"Input : {fragments}")
        print(f"Output: {outfile}")

        # ---- 1. Import fragments ----
        adata = snap.pp.import_fragments(
            fragment_file=fragments,
            chrom_sizes=snap.genome.mm39,
            sorted_by_barcode=False,
            n_jobs=8
        )
        adata.obs["sample"] = sample

        # ---- 2. Compute TSS enrichment per cell ----
        print("Computing TSS enrichment (snap.metrics.tsse)...")
        snap.metrics.tsse(adata, snap.genome.mm39)

        # ---- 3. Compute aggregate fragment size distribution ----
        print("Computing fragment size distribution...")
        size_counter = compute_fragment_size_distribution(fragments)

        # ---- 4. Compute average TSS profile ----
        print("Computing TSS profile...")
        tss_positions, tss_profile = compute_tss_profile(fragments, snap.genome.mm39)

        # ---- 5. Generate ONLY the unique QC plots (fragment size and TSS) ----
        print("Generating unique QC plots...")
        # Redundant histograms and scatter plots are removed because the filter script handles them (pre-filter and post-filter).
        plot_frag_size_dist_from_counter(size_counter, sample, "prefilter")
        plot_tss_profile(tss_positions, tss_profile, sample, "prefilter")

        # ---- 6. Save ----
        print("Saving...")
        adata.write(outfile)

        print(f"Finished {sample}")

        del adata, size_counter
        gc.collect()

    print("=" * 60)
    print("All samples imported successfully.")


if __name__ == "__main__":
    main()
