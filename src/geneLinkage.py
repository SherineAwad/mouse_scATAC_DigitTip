#!/usr/bin/env python3

import argparse
import os
import pickle

import numpy as np
import snapatac2 as snap


def convert_to_float32(adata, label):

    print()
    print(
        f"Loading {label} matrix into memory..."
    )

    adata_mem = adata.to_memory()

    print(
        f"{label} matrix loaded into memory."
    )

    print(
        f"{label} matrix type: {type(adata_mem.X)}"
    )

    print(
        f"{label} matrix dtype before conversion: "
        f"{adata_mem.X.dtype}"
    )

    adata_mem.X = adata_mem.X.astype(
        np.float32
    )

    print(
        f"{label} matrix dtype after conversion: "
        f"{adata_mem.X.dtype}"
    )

    return adata_mem


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Perform SnapATAC2 peak-to-gene linkage analysis "
            "using peak accessibility and gene activity matrices."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Peak-by-cell H5AD file."
    )

    parser.add_argument(
        "--gene_matrix",
        required=True,
        help="Cell-by-gene activity H5AD file."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output H5AD file."
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output files."
    )

    parser.add_argument(
        "--upstream",
        type=int,
        default=100000,
        help="Distance upstream of TSS. Default: 100000 bp."
    )

    parser.add_argument(
        "--downstream",
        type=int,
        default=100000,
        help="Distance downstream of TSS. Default: 100000 bp."
    )

    args = parser.parse_args()

    os.makedirs(
        "figures",
        exist_ok=True
    )

    # ------------------------------------------------------------
    # 1. Load peak matrix
    # ------------------------------------------------------------

    print()
    print("Loading peak matrix...")
    print(f"Input: {args.input}")

    peak_mat = snap.read(
        args.input
    )

    print(
        f"Peak matrix: "
        f"{peak_mat.n_obs} cells x "
        f"{peak_mat.n_vars} peaks"
    )

    # ------------------------------------------------------------
    # 2. Load gene activity matrix
    # ------------------------------------------------------------

    print()
    print("Loading gene activity matrix...")
    print(f"Gene matrix: {args.gene_matrix}")

    gene_mat = snap.read(
        args.gene_matrix
    )

    print(
        f"Gene matrix: "
        f"{gene_mat.n_obs} cells x "
        f"{gene_mat.n_vars} genes"
    )

    # ------------------------------------------------------------
    # 3. Check cell identities before conversion
    # ------------------------------------------------------------

    print()
    print("Checking cell identities...")

    peak_cells = list(
        peak_mat.obs_names
    )

    gene_cells = list(
        gene_mat.obs_names
    )

    peak_cell_set = set(
        peak_cells
    )

    gene_cell_set = set(
        gene_cells
    )

    if peak_cell_set != gene_cell_set:

        missing_from_gene = (
            peak_cell_set - gene_cell_set
        )

        missing_from_peak = (
            gene_cell_set - peak_cell_set
        )

        raise ValueError(
            "Peak matrix and gene matrix do not contain "
            "the same cells.\n"
            f"Cells in peak matrix but missing from gene matrix: "
            f"{len(missing_from_gene)}\n"
            f"Cells in gene matrix but missing from peak matrix: "
            f"{len(missing_from_peak)}"
        )

    if peak_cells != gene_cells:

        print(
            "Cell sets match but cell order differs."
        )

        print(
            "Reordering gene matrix to match peak matrix..."
        )

        gene_mat = gene_mat[
            peak_cells,
            :
        ].copy()

    print(
        "Cell identities and order match."
    )

    # ------------------------------------------------------------
    # 4. Read peak coordinates
    # ------------------------------------------------------------

    print()
    print("Reading peak coordinates...")

    peak_regions = list(
        peak_mat.var_names
    )

    if len(peak_regions) == 0:

        raise ValueError(
            "Peak matrix contains no peaks."
        )

    print(
        f"Candidate regulatory regions: "
        f"{len(peak_regions)}"
    )

    # ------------------------------------------------------------
    # 5. Build peak-to-gene candidate network
    # ------------------------------------------------------------

    print()
    print(
        "Building peak-to-gene candidate network..."
    )

    print(
        "Genome annotation: mm39"
    )

    print(
        f"Upstream distance: "
        f"{args.upstream} bp"
    )

    print(
        f"Downstream distance: "
        f"{args.downstream} bp"
    )

    network = snap.tl.init_network_from_annotation(
        regions=peak_regions,
        anno_file=snap.genome.mm39,
        upstream=args.upstream,
        downstream=args.downstream,
        id_type="gene_name",
        coding_gene_only=True
    )

    print(
        f"Network nodes: "
        f"{network.num_nodes()}"
    )

    print(
        f"Network edges: "
        f"{network.num_edges()}"
    )

    # ------------------------------------------------------------
    # 6. Convert backed CSR uint32 matrices to in-memory float32
    # ------------------------------------------------------------

    peak_mat_float = convert_to_float32(
        peak_mat,
        "Peak"
    )

    gene_mat_float = convert_to_float32(
        gene_mat,
        "Gene"
    )

    # ------------------------------------------------------------
    # 7. Calculate peak-gene correlation scores
    # ------------------------------------------------------------

    print()
    print(
        "Calculating peak-gene correlation scores..."
    )

    snap.tl.add_cor_scores(
        network,
        peak_mat=peak_mat_float,
        gene_mat=gene_mat_float,
        overwrite=True
    )

    print(
        "Correlation scores calculated."
    )

    # ------------------------------------------------------------
    # 8. Save network
    # ------------------------------------------------------------

    print()
    print(
        "Saving linkage network..."
    )

    network_file = (
        f"{args.prefix}_network.pkl"
    )

    with open(
        network_file,
        "wb"
    ) as f:

        pickle.dump(
            network,
            f
        )

    print(
        f"Saved network: {network_file}"
    )

    # ------------------------------------------------------------
    # 9. Save output H5AD
    # ------------------------------------------------------------

    print()
    print(
        f"Saving output H5AD: {args.output}"
    )

    peak_mat_float.write(
        args.output
    )

    print(
        f"Saved: {args.output}"
    )

    # ------------------------------------------------------------
    # 10. Done
    # ------------------------------------------------------------

    print()
    print(
        "Gene linkage analysis complete."
    )


if __name__ == "__main__":
    main()
