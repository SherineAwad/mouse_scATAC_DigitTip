#!/usr/bin/env python3

import argparse
import os
import pickle
import re

import pandas as pd
import matplotlib.pyplot as plt


def get_attribute(obj, key):
    """Get an attribute from a SnapATAC2/rustworkx node or edge."""
    if isinstance(obj, dict):
        return obj.get(key)
    return getattr(obj, key, None)


def extract_peak_name(node):
    """
    Extract genomic peak coordinate (chr:start-end).
    """
    chrom = get_attribute(node, "chrom") or get_attribute(node, "chr")
    start = get_attribute(node, "start")
    end = get_attribute(node, "end")

    if chrom is not None and start is not None and end is not None:
        return f"{chrom}:{start}-{end}"

    possible_keys = ["name", "id", "region", "peak"]
    values = []

    for key in possible_keys:
        value = get_attribute(node, key)
        if value is not None:
            values.append(str(value))

    pattern = r"(chr[^:\s]+:\d+-\d+)"

    for value in values:
        match = re.search(pattern, value)
        if match:
            return match.group(1)

    node_string = str(node)
    match = re.search(pattern, node_string)
    if match:
        return match.group(1)

    raise ValueError(f"Could not extract peak coordinate from node: {node}")


def extract_gene_name(node):
    """
    Extract gene name from node.
    """
    possible_keys = ["name", "id", "gene_name", "gene"]

    for key in possible_keys:
        value = get_attribute(node, key)
        if value is not None:
            value = str(value)
            if "|" in value:
                value = value.split("|")[0]
            return value

    node_string = str(node)
    patterns = [
        r"gene_name=['\"]([^'\"]+)['\"]",
        r"name=['\"]([^'\"]+)['\"]",
        r"id=['\"]([^'\"]+)['\"]"
    ]

    for pattern in patterns:
        match = re.search(pattern, node_string)
        if match:
            return match.group(1)

    raise ValueError(f"Could not extract gene name from node: {node}")


def get_node_type(node):
    """Get SnapATAC2 node type."""
    value = get_attribute(node, "type")
    if value is None:
        return None
    return str(value)


def get_correlation(edge_data):
    """
    Extract numeric correlation score from SnapATAC2 LinkData or Dict.
    """
    keys = ["cor_score", "correlation", "cor", "score"]
    
    for key in keys:
        value = get_attribute(edge_data, key)
        if value is not None:
            if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
                value = value[0]
            try:
                return float(value)
            except (ValueError, TypeError):
                continue

    return None


def extract_peak_gene_links(network):
    rows = []
    print("Extracting peak-gene correlations...")

    for source_idx, target_idx in network.edge_list():
        source_node = network[source_idx]
        target_node = network[target_idx]

        source_type = get_node_type(source_node)
        target_type = get_node_type(target_node)

        peak_node = None
        gene_node = None

        if source_type in ["region", "peak"] and target_type == "gene":
            peak_node = source_node
            gene_node = target_node
        elif target_type in ["region", "peak"] and source_type == "gene":
            peak_node = target_node
            gene_node = source_node
        else:
            continue

        edge_data = network.get_edge_data(source_idx, target_idx)
        correlation = get_correlation(edge_data)

        if correlation is None:
            continue

        peak = extract_peak_name(peak_node)
        gene = extract_gene_name(gene_node)

        rows.append(
            {
                "peak": peak,
                "gene": gene,
                "correlation": correlation
            }
        )

    if not rows:
        raise ValueError("No peak-gene correlations were found.")

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Plot peak-gene correlations from a SnapATAC2 network PKL."
    )
    parser.add_argument("--pkl", required=True, help="Peak-gene network PKL.")
    parser.add_argument("--prefix", required=True, help="Prefix for output files.")
    parser.add_argument(
        "--n",
        type=int,
        default=30,
        help="Number of strongest peak-gene pairs to plot. Default: 30."
    )

    args = parser.parse_args()

    if args.n < 1:
        raise ValueError("--n must be >= 1.")

    print(f"\nLoading network: {args.pkl}")
    with open(args.pkl, "rb") as f:
        network = pickle.load(f)

    print("Network loaded.")
    print(f"Nodes: {network.num_nodes()}")
    print(f"Edges: {network.num_edges()}")

    df = extract_peak_gene_links(network)
    print(f"Peak-gene pairs found: {len(df)}")

    df = df.drop_duplicates(subset=["peak", "gene"]).reset_index(drop=True)
    print(f"Unique peak-gene pairs: {len(df)}")

    df["abs_correlation"] = df["correlation"].abs()
    df = df.sort_values("abs_correlation", ascending=False).reset_index(drop=True)

    csv_output = f"{args.prefix}_peak_gene_correlations.csv"
    df.to_csv(csv_output, index=False)
    print(f"\nSaved CSV: {csv_output}")

    top_df = df.head(args.n).copy()
    top_df = top_df.sort_values("correlation").reset_index(drop=True)

    os.makedirs("figures", exist_ok=True)
    figure_output = os.path.join("figures", f"{args.prefix}_peak_gene_correlation.png")

    n_pairs = len(top_df)
    fig_height = max(6, n_pairs * 0.32 + 2)

    fig, ax = plt.subplots(figsize=(11, fig_height))

    y = list(range(n_pairs))

    # Clean, uniform marker sizes
    ax.scatter(top_df["correlation"], y, s=80, alpha=0.85, edgecolors="none")
    ax.axvline(0, linestyle="--", linewidth=1, color="gray")

    for position in y:
        ax.axhline(position, linewidth=0.4, alpha=0.15)

    # Neutral association notation (en-dash) instead of directional arrow
    labels = [f"{row['peak']} – {row['gene']}" for _, row in top_df.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)

    ax.set_xlabel("Peak-gene correlation")
    ax.set_ylabel("Peak – Gene Pair")
    ax.set_title(f"Top {n_pairs} peak-gene correlations")

    # Dynamic x-axis bounds fitting actual correlation range
    min_corr = top_df["correlation"].min()
    max_corr = top_df["correlation"].max()

    if min_corr > 0:
        x_min = min(0.0, min_corr - 0.05)
    else:
        x_min = min_corr - 0.05

    x_max = max_corr + 0.05

    ax.set_xlim(x_min, x_max)

    fig.tight_layout()
    fig.savefig(figure_output, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved figure: {figure_output}\n")
    print("Peak-gene correlation analysis complete.")
    print(f"Total links: {len(df)}")
    print(f"Links plotted: {n_pairs}")


if __name__ == "__main__":
    main()
