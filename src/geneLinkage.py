#!/usr/bin/env python3

import argparse
import gzip
import os
import pickle
import numpy as np
import snapatac2 as snap


def convert_to_float32(adata, label):
    print()
    print(f"Loading {label} matrix into memory...")
    adata_mem = adata.to_memory()

    # Use log1p layer if available to ensure correlation runs on normalized values
    if "log1p" in adata_mem.layers:
        print(f"Assigning 'log1p' layer to .X for {label} matrix...")
        adata_mem.X = adata_mem.layers["log1p"].copy()

    print(f"{label} matrix loaded into memory.")
    print(f"{label} matrix type: {type(adata_mem.X)}")
    print(f"{label} matrix dtype before conversion: {adata_mem.X.dtype}")

    if adata_mem.X.dtype != np.float32:
        adata_mem.X = adata_mem.X.astype(np.float32, copy=False)

    print(f"{label} matrix dtype after conversion: {adata_mem.X.dtype}")
    return adata_mem


def exclude_promoter_peaks(peak_regions, exclude_bp=2000):
    """
    Excludes peaks falling within +/- exclude_bp of protein-coding gene TSS coordinates.
    """
    print(f"Excluding peaks within +/- {exclude_bp} bp of protein-coding TSS...")

    # Extract TSS positions from mm39 gene annotation for protein-coding genes only
    anno_path = str(snap.genome.mm39.annotation)
    tss_dict = {}

    open_fn = gzip.open if anno_path.endswith(".gz") else open
    with open_fn(anno_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9 and parts[2] == "gene":
                attributes = parts[8]
                # Match protein-coding genes to stay consistent with coding_gene_only=True
                if "protein_coding" in attributes:
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    tss = start if strand in ["+", "1"] else end
                    if chrom not in tss_dict:
                        tss_dict[chrom] = []
                    tss_dict[chrom].append(tss)

    # Convert TSS lists to sorted numpy arrays for fast binary search
    for chrom in tss_dict:
        tss_dict[chrom] = np.array(sorted(tss_dict[chrom]))

    filtered_peaks = []
    excluded_count = 0

    for peak in peak_regions:
        try:
            chrom, coords = peak.split(":")
            p_start, p_end = map(int, coords.split("-"))
            p_mid = (p_start + p_end) // 2

            if chrom in tss_dict:
                idx = np.searchsorted(tss_dict[chrom], p_mid)
                distances = []
                if idx < len(tss_dict[chrom]):
                    distances.append(abs(tss_dict[chrom][idx] - p_mid))
                if idx > 0:
                    distances.append(abs(tss_dict[chrom][idx - 1] - p_mid))

                if min(distances) <= exclude_bp:
                    excluded_count += 1
                    continue

            filtered_peaks.append(peak)
        except Exception as e:
            raise ValueError(f"Failed to parse peak string '{peak}': {e}")

    print(f"Excluded {excluded_count} promoter-proximal peaks.")
    print(f"Remaining candidate distal peaks: {len(filtered_peaks)}")
    return filtered_peaks


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Perform SnapATAC2 peak-to-gene linkage analysis "
            "using peak accessibility and gene activity matrices."
        )
    )

    parser.add_argument(
        "--input", required=True, help="Peak-by-cell H5AD file."
    )

    parser.add_argument(
        "--gene_matrix", required=True, help="Cell-by-gene activity H5AD file."
    )

    parser.add_argument(
        "--prefix",
        required=False,
        default=None,
        help="Prefix for output network file. If omitted, derives from --output.",
    )

    parser.add_argument("--output", required=True, help="Output H5AD file.")

    parser.add_argument(
        "--upstream",
        type=int,
        default=100000,
        help="Distance upstream of TSS. Default: 100000 bp.",
    )

    parser.add_argument(
        "--downstream",
        type=int,
        default=100000,
        help="Distance downstream of TSS. Default: 100000 bp.",
    )

    parser.add_argument(
        "--exclude-promoter-bp",
        type=int,
        default=2000,
        help="Exclusion window around TSS in bp. Default: 2000 (+/- 2kb).",
    )

    args = parser.parse_args()

    # Determine prefix
    prefix = args.prefix if args.prefix else os.path.splitext(args.output)[0]

    # ------------------------------------------------------------
    # 1. Load peak matrix
    # ------------------------------------------------------------

    print()
    print("Loading peak matrix...")
    print(f"Input: {args.input}")

    peak_mat = snap.read(args.input, backed="r")

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

    gene_mat = snap.read(args.gene_matrix, backed="r")

    print(
        f"Gene matrix: "
        f"{gene_mat.n_obs} cells x "
        f"{gene_mat.n_vars} genes"
    )

    # ------------------------------------------------------------
    # 3. Check cell identities and subset gene matrix
    # ------------------------------------------------------------

    print()
    print("Checking cell identities...")

    # Load into memory and resolve duplicate barcodes
    peak_mat_mem = peak_mat.to_memory()
    gene_mat_mem = gene_mat.to_memory()

    peak_mat_mem.obs_names_make_unique()
    gene_mat_mem.obs_names_make_unique()

    peak_cells = list(peak_mat_mem.obs_names)
    gene_obs_names = list(gene_mat_mem.obs_names)

    gene_cell_dict = {name: idx for idx, name in enumerate(gene_obs_names)}

    missing_cells = [cell for cell in peak_cells if cell not in gene_cell_dict]
    if len(missing_cells) > 0:
        raise ValueError(
            "Peak matrix contains cells missing from gene matrix.\n"
            f"Cells missing: {len(missing_cells)}"
        )

    print(
        f"Subsetting gene matrix from {len(gene_obs_names)} cells "
        f"to match {len(peak_cells)} peak cells..."
    )

    # Slice by positional integer indices
    target_indices = [gene_cell_dict[cell] for cell in peak_cells]
    gene_mat = gene_mat_mem[target_indices, :].copy()
    peak_mat = peak_mat_mem

    print("Cell identities successfully matched and ordered.")

    # ------------------------------------------------------------
    # 4. Read peak coordinates & Exclude +/- 2kb Promoter Regions
    # ------------------------------------------------------------

    print()
    print("Reading peak coordinates...")

    peak_regions = list(peak_mat.var_names)

    if len(peak_regions) == 0:
        raise ValueError("Peak matrix contains no peaks.")

    print(f"Total input regions: {len(peak_regions)}")

    distal_peak_regions = exclude_promoter_peaks(
        peak_regions, exclude_bp=args.exclude_promoter_bp
    )

    # ------------------------------------------------------------
    # 5. Build peak-to-gene candidate network
    # ------------------------------------------------------------

    print()
    print("Building peak-to-gene candidate network...")
    print("Genome annotation: mm39")
    print(f"Upstream distance: {args.upstream} bp")
    print(f"Downstream distance: {args.downstream} bp")

    network = snap.tl.init_network_from_annotation(
        regions=distal_peak_regions,
        anno_file=snap.genome.mm39,
        upstream=args.upstream,
        downstream=args.downstream,
        id_type="gene_name",
        coding_gene_only=True,
    )

    print(f"Network nodes: {network.num_nodes()}")
    print(f"Network edges: {network.num_edges()}")

    # ------------------------------------------------------------
    # 6. Convert backed CSR uint32 matrices to in-memory float32
    # ------------------------------------------------------------

    peak_mat_float = convert_to_float32(peak_mat, "Peak")
    gene_mat_float = convert_to_float32(gene_mat, "Gene")

    # ------------------------------------------------------------
    # 7. Calculate peak-gene correlation scores
    # ------------------------------------------------------------

    print()
    print("Calculating peak-gene correlation scores...")

    snap.tl.add_cor_scores(
        network,
        peak_mat=peak_mat_float,
        gene_mat=gene_mat_float,
        overwrite=True,
    )

    print("Correlation scores calculated.")

    # ------------------------------------------------------------
    # 8. Save peak-gene network as PKL
    # ------------------------------------------------------------

    network_output = f"{prefix}_network.pkl"

    print()
    print(f"Saving peak-gene network: {network_output}")

    with open(network_output, "wb") as f:
        pickle.dump(network, f)

    print(f"Saved network: {network_output}")

    # ------------------------------------------------------------
    # 9. Save output H5AD
    # ------------------------------------------------------------

    print()
    print(f"Saving output H5AD: {args.output}")

    peak_mat_float.write(args.output)

    print(f"Saved: {args.output}")
    print("Gene linkage analysis complete.")


if __name__ == "__main__":
    main()
