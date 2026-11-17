#!/usr/bin/env python3

import argparse
import gc
import re
from pathlib import Path

import anndata as ad
import circe as ci
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import scipy.sparse as sp


# ----------------------------------------------------------------------
# Direct h5py readers for the X matrix in an .h5ad file
# ----------------------------------------------------------------------

def _load_cells_chunk(h5_file, cell_indices, n_vars):
    g = h5_file["X"]

    if not isinstance(g, h5py.Group):
        return sp.csr_matrix(g[...][cell_indices, :])

    encoding = g.attrs.get("encoding-type", "csr_matrix")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    enc = encoding.lower()

    cell_indices = np.asarray(cell_indices, dtype=np.int64)
    n_rows = len(cell_indices)

    if enc in ("csr_matrix", "csr"):
        indptr = g["indptr"][:]
        data_ds, indices_ds = g["data"], g["indices"]
        sort_order = np.argsort(cell_indices, kind="stable")
        sorted_cells = cell_indices[sort_order]
        starts, ends = indptr[sorted_cells], indptr[sorted_cells + 1]
        sizes = ends - starts
        dtype = data_ds.dtype
        if sizes.sum() == 0:
            return sp.csr_matrix((n_rows, n_vars), dtype=dtype)
        positions = np.concatenate(
            [np.arange(s, e, dtype=np.int64) for s, e in zip(starts, ends)]
        )
        all_data = data_ds[positions]
        all_indices = indices_ds[positions]
        indptr2 = np.zeros(n_rows + 1, dtype=np.int64)
        np.cumsum(sizes, out=indptr2[1:])
        mat = sp.csr_matrix((all_data, all_indices, indptr2),
                            shape=(n_rows, n_vars))
        inv = np.empty_like(sort_order)
        inv[sort_order] = np.arange(n_rows)
        return mat[inv]

    if enc in ("csc_matrix", "csc"):
        shape = g.attrs["shape"]
        csc = sp.csc_matrix(
            (g["data"][:], g["indices"][:], g["indptr"][:]),
            shape=(int(shape[0]), int(shape[1])),
        )
        return csc.tocsr()[cell_indices, :]

    raise ValueError(f"Unsupported X encoding-type: {encoding}")


def load_cells_from_h5ad(h5_path, cell_indices, n_vars, batch_size=20000):
    cell_indices = np.asarray(cell_indices, dtype=np.int64)
    with h5py.File(h5_path, "r") as f:
        if len(cell_indices) <= batch_size:
            return _load_cells_chunk(f, cell_indices, n_vars)
        mats = []
        for i in range(0, len(cell_indices), batch_size):
            mats.append(_load_cells_chunk(f, cell_indices[i:i + batch_size], n_vars))
            gc.collect()
        return sp.vstack(mats, format="csr")


# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Condition-specific CIRCE co-accessibility analysis"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--genome", default="mouse",
                        choices=["mouse", "human", "drosophila"])

    # Memory-control knobs
    parser.add_argument("--max-regions", type=int, default=200_000,
                        help="Cap regions per condition (keep top by "
                             "accessibility). Default 200000.")
    parser.add_argument("--max-cells", type=int, default=20_000,
                        help="Subsample cells per condition before metacells "
                             "if greater than this. Set <=0 to disable.")
    parser.add_argument("--min-cells", type=int, default=30,
                        help="Minimum number of cells a region must be "
                             "present in (default: 30).")
    parser.add_argument("--min-cell-frac", type=float, default=0.01,
                        help="Minimum fraction of cells a region must be "
                             "present in (default: 0.01). Set 0 to disable.")
    parser.add_argument("--min-condition-cells", type=int, default=200,
                        help="Skip conditions with fewer cells than this after "
                             "subsampling (default: 200). Set 0 to disable.")

    args = parser.parse_args()

    Path("figures").mkdir(exist_ok=True)

    print("=" * 65, flush=True)
    print("CIRCE: Condition-specific co-accessibility analysis", flush=True)
    print("=" * 65, flush=True)
    print(f"\nInput:  {args.input}", flush=True)
    print(f"Genome: {args.genome}", flush=True)
    print(f"max-regions: {args.max_regions:,}", flush=True)
    print(f"max-cells:   {args.max_cells if args.max_cells > 0 else 'disabled'}",
          flush=True)

    # --------------------------------------------------
    # 1. Metadata only
    # --------------------------------------------------
    print("\nOpening AnnData (metadata only)...", flush=True)
    adata = ad.read_h5ad(args.input, backed="r")
    n_vars = adata.n_vars
    print(f"Cells: {adata.n_obs:,}", flush=True)
    print(f"Regions: {n_vars:,}", flush=True)

    if "sample" not in adata.obs.columns:
        raise ValueError("Missing 'sample' column in adata.obs")
    missing = [c for c in ("chromosome", "start", "end")
               if c not in adata.var.columns]
    if missing:
        raise ValueError(f"Missing genomic coordinate columns: {missing}")

    obs_df = adata.obs.copy()
    var_df = adata.var.copy()
    adata.file.close()
    del adata
    gc.collect()

    # --------------------------------------------------
    # 2. Conditions
    # --------------------------------------------------
    conditions = obs_df["sample"].dropna().unique().tolist()
    print(f"\nConditions found: {len(conditions)}", flush=True)
    for c in conditions:
        n = (obs_df["sample"] == c).sum()
        print(f"  {c}: {n:,} cells", flush=True)

    # --------------------------------------------------
    # 3. Per-condition processing
    # --------------------------------------------------
    for condition in conditions:
        print("\n" + "=" * 65, flush=True)
        print(f"Processing condition: {condition}", flush=True)
        print("=" * 65, flush=True)

        safe_condition = re.sub(r"[^A-Za-z0-9._-]+", "_", str(condition)).strip("_")

        # -------- 3.1 select (and maybe subsample) cells
        cell_mask = obs_df["sample"].to_numpy() == condition
        cell_indices = np.flatnonzero(cell_mask)
        n_cond = len(cell_indices)
        if n_cond == 0:
            print("  Skipping: no cells.", flush=True)
            continue

        if args.max_cells > 0 and n_cond > args.max_cells:
            print(
                f"  Subsampling cells for metacells: {n_cond:,} -> "
                f"{args.max_cells:,} (fixed seed=0)",
                flush=True,
            )
            rng = np.random.default_rng(0)
            pick = np.sort(rng.choice(n_cond, args.max_cells, replace=False))
            cell_indices = cell_indices[pick]
            n_cond = args.max_cells

        if args.min_condition_cells > 0 and n_cond < args.min_condition_cells:
            print(f"  Skipping: only {n_cond} cells "
                  f"(< --min-condition-cells={args.min_condition_cells})", flush=True)
            continue

        print(f"\nLoading matrix for {n_cond:,} cells × {n_vars:,} regions...",
              flush=True)
        mat_cells = load_cells_from_h5ad(args.input, cell_indices, n_vars)
        print(f"  Loaded: nnz={mat_cells.nnz:,}", flush=True)

        # -------- 3.2 filter + cap regions
        print("\nFiltering regions...", flush=True)
        min_cells_thresh = max(
            args.min_cells,
            int(round(args.min_cell_frac * n_cond)),
        )
        print(f"  min cells per region (per-condition): {min_cells_thresh}",
              flush=True)

        chunk_size = 500_000
        valid_vars_list, valid_sum_list = [], []
        for s in range(0, n_vars, chunk_size):
            e = min(s + chunk_size, n_vars)
            ch = mat_cells[:, s:e]
            counts = np.asarray(ch.getnnz(axis=0)).flatten()
            sums   = np.asarray(ch.sum(axis=0)).flatten()
            mask   = counts >= min_cells_thresh
            valid_vars_list.append(np.flatnonzero(mask) + s)
            valid_sum_list.append(sums[mask])
            del ch
            gc.collect()

        var_indices = (np.concatenate(valid_vars_list)
                       if valid_vars_list else np.array([], np.int64))
        var_sums = (np.concatenate(valid_sum_list)
                    if valid_sum_list else np.array([], np.float64))
        assert len(var_indices) == len(var_sums), "index/sum desync"
        print(f"  Passed filter: {len(var_indices):,}", flush=True)

        if len(var_indices) == 0:
            print("  Skipping: no regions retained.", flush=True)
            del mat_cells
            gc.collect()
            continue

        if len(var_indices) > args.max_regions:
            print(
                f"  Capping to top {args.max_regions:,} by total accessibility",
                flush=True,
            )
            top = np.argsort(var_sums)[::-1][:args.max_regions]
            var_indices = np.sort(var_indices[top])

        print(f"  Regions retained: {len(var_indices):,}", flush=True)
        del var_sums, valid_vars_list, valid_sum_list
        gc.collect()

        # -------- 3.3 build condition AnnData
        condition_mat = mat_cells[:, var_indices].tocsr()
        del mat_cells
        gc.collect()

        condition_adata = ad.AnnData(
            X=condition_mat,
            obs=obs_df.iloc[cell_indices].copy(),
            var=var_df.iloc[var_indices].copy(),
        )
        print(f"  condition_adata: {condition_adata.shape}, "
              f"nnz={condition_adata.X.nnz:,}", flush=True)

        if not sparse.isspmatrix_csr(condition_adata.X):
            condition_adata.X = condition_adata.X.tocsr()

        # -------- 3.4 metacells
        print("\nComputing metacells...", flush=True)
        metacells = ci.metacells.compute_metacells(condition_adata)
        print(f"  Metacells: {metacells.n_obs:,}", flush=True)
        if metacells.n_obs < 20:
            print(f"  WARNING: only {metacells.n_obs} metacells — "
                  f"co-accessibility estimates will be noisy.", flush=True)
        print(f"  Regions:   {metacells.n_vars:,}", flush=True)
        print(f"  metacells.X type: {type(metacells.X).__name__}, "
              f"shape={metacells.X.shape}", flush=True)

        # -------- 3.5 release original
        print("\nReleasing original condition matrix...", flush=True)
        del condition_adata, cell_indices, cell_mask, var_indices, condition_mat
        gc.collect()

        # -------- 3.6 CSC (sparse/dense safe)
        print("\nConverting metacell matrix to CSC...", flush=True)
        if sparse.issparse(metacells.X):
            metacells.X = metacells.X.tocsc()
        else:
            metacells.X = sp.csc_matrix(np.asarray(metacells.X))
        gc.collect()

        # -------- 3.7 network
        print("\nComputing co-accessibility network...", flush=True)
        ci.compute_atac_network(metacells, organism=args.genome,
                                njobs=1, verbose=1)
        print("\nCo-accessibility network completed.", flush=True)

        # -------- 3.8 links
        print("\nExtracting co-accessibility links...", flush=True)
        links = ci.extract_atac_links(metacells)
        links_file = f"{args.prefix}_{safe_condition}_coaccessibility_links.csv"
        links.to_csv(links_file, index=False)
        print(f"  Links: {len(links):,}", flush=True)
        print(f"  Saved: {links_file}", flush=True)

        # -------- 3.9 CCANs
        print("\nIdentifying CCANs...", flush=True)
        if not links.empty:
            ccans = ci.find_ccans(links, seed=0)
        else:
            ccans = pd.DataFrame(columns=["Peak", "CCAN"])
        ccans_file = f"{args.prefix}_{safe_condition}_ccans.csv"
        ccans.to_csv(ccans_file, index=False)
        print(f"  CCAN assignments: {len(ccans):,}", flush=True)
        print(f"  Saved: {ccans_file}", flush=True)

        del metacells, links, ccans
        gc.collect()
        print(f"\nFinished condition: {condition}", flush=True)

    print("\n" + "=" * 65, flush=True)
    print("CIRCE ANALYSIS COMPLETED", flush=True)
    print("=" * 65, flush=True)
    print(f"Conditions processed: {len(conditions)}", flush=True)
    print(f"Genome: {args.genome}", flush=True)
    print(f"Output prefix: {args.prefix}", flush=True)


if __name__ == "__main__":
    main()
