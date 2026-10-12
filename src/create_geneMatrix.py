
#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd
import scanpy as sc
import snapatac2 as snap


def parse_gff3_tss(gff3_path):
    genes = []
    open_fn = gzip.open if str(gff3_path).endswith(".gz") else open

    with open_fn(gff3_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            if parts[2] == "gene":
                chromosome = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                gene_name = None
                for attr in attributes.split(";"):
                    if attr.startswith("gene_name=") or attr.startswith("Name="):
                        gene_name = attr.split("=", 1)[1]
                        break

                if gene_name:
                    tss = start if strand in ["+", "1"] else end

                    genes.append({
                        "gene_name": gene_name,
                        "chromosome": chromosome,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "tss": tss,
                    })

    df = pd.DataFrame(genes)

    if df.empty:
        raise ValueError(f"No gene annotations found in {gff3_path}")

    return df.drop_duplicates(subset=["gene_name"]).set_index("gene_name")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create a normalized, log1p-transformed gene-activity "
            "matrix from a SnapATAC2 ATAC object."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input SnapATAC2 h5ad file",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output gene-matrix h5ad file",
    )

    args = parser.parse_args()

    print("Reading input object...")
    adata = snap.read(args.input)

    print(f"Loaded: {adata.n_obs} cells")

    print("Creating gene-activity matrix using mm39 annotation...")

    gene_matrix = snap.pp.make_gene_matrix(
        adata,
        gene_anno=snap.genome.mm39,
        inplace=False,
    )

    gene_matrix = (
        gene_matrix.to_memory()
        if hasattr(gene_matrix, "to_memory")
        else gene_matrix
    )

    print("Parsing mm39 TSS annotations...")

    anno_df = parse_gff3_tss(snap.genome.mm39.annotation)
    matched_var = anno_df.reindex(gene_matrix.var_names)

    for col in ["chromosome", "start", "end", "strand", "tss"]:
        gene_matrix.var[col] = matched_var[col].values

    print("Normalizing gene activity and applying log1p...")

    # Preserve the original gene-activity counts.
    gene_matrix.layers["counts"] = gene_matrix.X.copy()

    # Normalize each cell to 10,000 total gene-activity counts,
    # then apply log1p.
    sc.pp.normalize_total(gene_matrix, target_sum=1e4)
    sc.pp.log1p(gene_matrix)

    # Save the log-transformed matrix explicitly.
    gene_matrix.layers["log1p"] = gene_matrix.X.copy()

    print(
        f"Gene matrix: {gene_matrix.n_obs} cells x "
        f"{gene_matrix.n_vars} genes"
    )

    print("Saved layers: counts, log1p")
    print("Saving gene matrix...")

    gene_matrix.write(args.output)

    print(f"Saved: {args.output}")
    print("Done.")


if __name__ == "__main__":
    main()
