#!/usr/bin/env python3

import argparse
import anndata as ad

def main():
    parser = argparse.ArgumentParser(
        description="Rename samples in an AnnData object."
    )
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--output", required=True, help="Output h5ad file with renamed samples")
    args = parser.parse_args()

    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    # Check if 'sample' column exists
    if "sample" not in adata.obs.columns:
        raise ValueError("Column 'sample' not found in adata.obs")

    # Rename mapping
    rename_map = {
        "SINAA6": "control",
        "SINAA8": "3xAmp",
        "SINAB6": "1xAmp",
        "SINAH5": "5xAmp"
    }

    # Apply renaming
    adata.obs["sample"] = adata.obs["sample"].replace(rename_map)

    print(f"Sample names updated. Unique samples: {adata.obs['sample'].unique().tolist()}")

    print(f"Saving to {args.output}...")
    adata.write(args.output, compression="gzip")
    print("Done.")

if __name__ == "__main__":
    main()
