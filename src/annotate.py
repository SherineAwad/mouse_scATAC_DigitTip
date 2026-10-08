#!/usr/bin/env python3

import argparse
import os

import pandas as pd
import snapatac2 as snap
import matplotlib.pyplot as plt


def read_annotations(annotation_file):

    annotations = {}

    with open(annotation_file, "r") as f:
        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            cluster, annotation = line.split(",", 1)

            annotations[str(cluster).strip()] = annotation.strip()

    return annotations


def main():

    parser = argparse.ArgumentParser(
        description="Annotate Leiden clusters in a gene activity H5AD."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input geneMatrix H5AD file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output annotated H5AD file"
    )

    parser.add_argument(
        "--annotations",
        required=True,
        help="Cluster annotation file: cluster,annotation"
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output figures"
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # 1. Read gene matrix
    # ------------------------------------------------------------

    print("Reading gene matrix...")

    adata = snap.read(args.input)
    adata = adata.to_memory()

    print(
        f"Loaded: "
        f"{adata.n_obs} cells x "
        f"{adata.n_vars} genes"
    )

    # ------------------------------------------------------------
    # 2. Check Leiden clusters
    # ------------------------------------------------------------

    if "leiden" not in adata.obs.columns:
        raise ValueError(
            "Input file does not contain adata.obs['leiden']"
        )

    print(
        f"Leiden clusters: "
        f"{adata.obs['leiden'].nunique()}"
    )

    # ------------------------------------------------------------
    # 3. Read cluster annotations
    # ------------------------------------------------------------

    print(
        f"Reading annotations from {args.annotations}..."
    )

    annotations = read_annotations(args.annotations)

    print(
        f"Annotations provided for "
        f"{len(annotations)} clusters"
    )

    # ------------------------------------------------------------
    # 4. Check that all Leiden clusters are annotated
    # ------------------------------------------------------------

    clusters = set(
        adata.obs["leiden"].astype(str).unique()
    )

    missing = clusters - set(annotations.keys())
    extra = set(annotations.keys()) - clusters

    if missing:
        raise ValueError(
            "Missing annotations for clusters: "
            + ", ".join(sorted(missing, key=int))
        )

    if extra:
        print(
            "Warning: annotations provided for clusters "
            "not present in the H5AD: "
            + ", ".join(sorted(extra, key=int))
        )

    # ------------------------------------------------------------
    # 5. Add annotation to obs
    # ------------------------------------------------------------

    adata.obs["annotation"] = (
        adata.obs["leiden"]
        .astype(str)
        .map(annotations)
    )

    print("\nCluster annotations:")

    summary = (
        adata.obs[["leiden", "annotation"]]
        .drop_duplicates()
        .sort_values("leiden", key=lambda x: x.astype(int))
    )

    print(summary.to_string(index=False))

    # ------------------------------------------------------------
    # 6. Save annotated H5AD
    # ------------------------------------------------------------

    print(
        f"\nSaving annotated H5AD to {args.output}..."
    )

    adata.write(args.output)

    print(f"Saved: {args.output}")

    # ------------------------------------------------------------
    # 7. Plot annotated UMAP
    # ------------------------------------------------------------

    os.makedirs("figures", exist_ok=True)

    if "X_umap" in adata.obsm:

        print("Creating annotated UMAP...")

        plt.figure(figsize=(10, 8))

        categories = adata.obs["annotation"].astype("category")
        coords = adata.obsm["X_umap"]

        for annotation in categories.cat.categories:

            mask = categories == annotation

            plt.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=2,
                label=annotation
            )

        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title("Cell type annotation")
        plt.legend(
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            markerscale=4
        )

        plt.tight_layout()

        figure_file = (
            f"figures/{args.prefix}_umap.png"
        )

        plt.savefig(
            figure_file,
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

        print(f"Saved: {figure_file}")

    else:

        print(
            "No X_umap found in input. "
            "Skipping UMAP plot."
        )

    print("\nDone.")

    # ------------------------------------------------------------
    # Annotated UMAP with cell-type labels
    # ------------------------------------------------------------

    print("Creating annotated UMAP with cell-type labels...")

    os.makedirs("figures", exist_ok=True)

    if "X_umap" not in adata.obsm:
        raise ValueError(
        "Input file does not contain X_umap."
    )

    coords = adata.obsm["X_umap"]

    plt.figure(figsize=(12, 10))

    annotations = adata.obs["annotation"].astype("category")

    for annotation in annotations.cat.categories:

        mask = annotations == annotation

        plt.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=2,
            alpha=0.7,
            label=annotation
        )

        # Centre of the cells belonging to this annotation
        x = coords[mask, 0].mean()
        y = coords[mask, 1].mean()

        plt.text(
            x,
            y,
            annotation,
            fontsize=10,
            fontweight="bold",
            ha="center",
            va="center",
            bbox=dict(
                facecolor="white",
                alpha=0.7,
                edgecolor="none"
             )
         )

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Cell-type annotation")

    plt.tight_layout()

    figure_file = (
    f"figures/{args.prefix}_umapON.png"
    )

    plt.savefig(
        figure_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Saved: {figure_file}")

if __name__ == "__main__":
    main()
