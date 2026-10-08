import argparse
import numpy as np
import pychromvar as pc
import snapatac2 as snap
from pyjaspar import jaspardb


class MotifWrapper:

    def __init__(self, motif):
        self._motif = motif

    @property
    def matrix_id(self):
        return getattr(
            self._motif, "id", getattr(self._motif, "name", "motif")
        )

    @property
    def counts(self):
        freqs = getattr(
            self._motif,
            "frequencies",
            getattr(self._motif, "matrix", None)
        )

        if freqs is None:
            raise AttributeError(
                "Motif object lacks frequencies or matrix data."
            )

        mat = np.array(freqs)

        if mat.shape[0] == 4 and mat.shape[1] != 4:
            mat = mat.T

        counts_mat = (mat * 1000).astype(int)

        return {
            "A": counts_mat[:, 0].tolist(),
            "C": counts_mat[:, 1].tolist(),
            "G": counts_mat[:, 2].tolist(),
            "T": counts_mat[:, 3].tolist(),
        }

    def __getattr__(self, attr):
        return getattr(self._motif, attr)


def main():

    parser = argparse.ArgumentParser(
        description="Run pychromVAR on AnnData"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to input h5ad file"
    )

    parser.add_argument(
        "--output",
        default="chromVAR.h5ad",
        help="Path to output file"
    )

    parser.add_argument(
        "--prefix",
        default="chromVAR",
        help="Prefix for output"
    )

    parser.add_argument(
        "--genome",
        default="mm39",
        help="Genome build"
    )

    parser.add_argument(
        "--motifs",
        default="JASPAR2020",
        help="JASPAR release"
    )

    parser.add_argument(
        "--n-peaks",
        type=int,
        default=50000,
        help="Maximum number of peaks to use"
    )

    parser.add_argument(
        "--n-cpu",
        type=int,
        default=4,
        help="Number of CPUs/jobs to use"
    )

    args = parser.parse_args()

    print(f"Loading AnnData from {args.input}...")
    adata = snap.read(args.input, backed=None)

    if adata.n_vars > args.n_peaks:
        print(
            f"Subsetting peaks from {adata.n_vars} "
            f"to {args.n_peaks}..."
        )

        top_peaks = np.argsort(
            np.array(adata.X.sum(axis=0)).flatten()
        )[-args.n_peaks:]

        adata = adata[:, top_peaks].copy()

    print(
        f"Fetching motifs ({args.motifs}) "
        f"and applying MotifWrapper..."
    )

    jdb_obj = jaspardb(release=args.motifs)

    raw_motifs = jdb_obj.fetch_motifs(
        collection="CORE",
        tax_group=["vertebrates"]
    )

    motifs = [MotifWrapper(m) for m in raw_motifs]

    print(f"Loaded {len(motifs)} motifs.")

    print(
        f"Adding peak sequences using SnapATAC2 {args.genome}..."
    )

    genome_file = snap.genome.mm39.fasta

    pc.add_peak_seq(
        adata,
        genome_file=genome_file,
        delimiter=r"[-:]"
    )

    print("Calculating GC bias...")
    pc.add_gc_bias(adata)

    print("Finding background regions...")
    pc.get_bg_peaks(
        adata,
        n_jobs=args.n_cpu
    )

    print("Matching motifs...")
    pc.match_motif(
        adata,
        motifs
    )

    print("Converting count matrix to float32...")
    adata.X = adata.X.astype(np.float32)

    print("Computing chromVAR deviations...")
    dev = pc.compute_deviations(
        adata,
        n_jobs=args.n_cpu
    )

    print(f"Saving output matrix to {args.output}...")
    dev.write_h5ad(args.output)

    print("Done!")


if __name__ == "__main__":
    main()
