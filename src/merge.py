#!/usr/bin/env python3

import argparse
import os
import anndata as ad
import snapatac2 as snap


def read_samples(sample_file):
    with open(sample_file) as f:
        return [line.strip() for line in f if line.strip()]


def main():

    parser = argparse.ArgumentParser(
        description="Merge snapATAC2 h5ad sample objects."
    )

    parser.add_argument(
        "--samples",
        required=True,
        help="Text file containing one sample name per line."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output merged h5ad filename."
    )

    args = parser.parse_args()

    samples = read_samples(args.samples)

    adatas = []

    for sample in samples:

        infile = f"{sample}_filtered.h5ad"

        if not os.path.exists(infile):
            raise FileNotFoundError(
                f"Cannot find {infile}"
            )

        print(f"Loading {infile}")

        adata = ad.read_h5ad(infile)

        print(
            f"{sample}: {adata.shape[0]} cells, {adata.shape[1]} features"
        )

        adatas.append(adata)

    print(f"Loaded {len(adatas)} samples")

    print("Merging samples...")

    snap.concat(
        adatas,
        file=args.output
    )

    print(f"Saved merged object: {args.output}")


if __name__ == "__main__":
    main()
