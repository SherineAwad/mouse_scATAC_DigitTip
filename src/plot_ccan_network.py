import argparse
import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    parser = argparse.ArgumentParser(description="Plot CCAN co-accessibility network graph.")
    parser.add_argument("--links", required=True)
    parser.add_argument("--ccans", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cutoff", type=float, default=None)
    parser.add_argument("--n", type=int, default=None, help="Plot only top N CCAN modules.")
    parser.add_argument("--nlinks", type=int, default=5000, help="Maximum number of edges to plot.")
    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)

    print(f"Loading CCANs from {args.ccans}...", flush=True)
    ccans_df = pd.read_csv(args.ccans)

    print(f"Loading links from {args.links}...", flush=True)
    links_df = pd.read_csv(args.links)
    print(f"Total raw edges loaded: {len(links_df):,}", flush=True)
    
    if args.cutoff is not None:
        links_df = links_df[links_df['score'] >= args.cutoff]
        print(f"Filtered to score >= {args.cutoff}. Remaining edges: {len(links_df):,}", flush=True)

    # --- NEW GLOBAL SUMMARY BLOCK ---
    global_peaks = set(links_df['Peak1']).union(set(links_df['Peak2']))
    print("\n" + "="*50, flush=True)
    print(f"GLOBAL SUMMARY FOR: {args.prefix}", flush=True)
    print(f"Total edges passing cutoff: {len(links_df):,}", flush=True)
    print(f"Total unique peaks involved: {len(global_peaks):,}", flush=True)
    print("="*50 + "\n", flush=True)
    # --------------------------------

    if args.n is not None:
        print(f"Selecting top {args.n} modules by peak count...", flush=True)
        top_ccans = ccans_df['CCAN'].value_counts().nlargest(args.n).index.tolist()
        ccans_df = ccans_df[ccans_df['CCAN'].isin(top_ccans)]
        print(f"Retained CCANs: {top_ccans}", flush=True)

    valid_peaks = set(ccans_df['Peak'])
    links_df = links_df[links_df['Peak1'].isin(valid_peaks) & links_df['Peak2'].isin(valid_peaks)]
    
    # USING THE NEW nlinks PARAMETER HERE
    if len(links_df) > args.nlinks:
        links_df = links_df.sort_values(by='score', ascending=False).head(args.nlinks)
        print(f"Hard limit applied: Plotting top {args.nlinks} edges to prevent timeout.", flush=True)

    print(f"Final edges to plot: {len(links_df):,}", flush=True)

    node_to_ccan = dict(zip(ccans_df['Peak'], ccans_df['CCAN']))
    unique_ccans = ccans_df['CCAN'].unique()
    
    palette = sns.color_palette("tab20", n_colors=len(unique_ccans))
    ccan_color_map = {ccan: palette[i] for i, ccan in enumerate(unique_ccans)}

    G = nx.Graph()

    for peak, ccan in node_to_ccan.items():
        G.add_node(peak, ccan=ccan)

    for _, row in links_df.iterrows():
        G.add_edge(row['Peak1'], row['Peak2'], weight=row['score'])

    G.remove_nodes_from(list(nx.isolates(G)))
    print(f"Network built: {G.number_of_nodes()} connected nodes, {G.number_of_edges()} edges.", flush=True)

    if G.number_of_nodes() == 0:
        print("Error: No connected nodes left after filtering.", flush=True)
        return

    print("Calculating layout...", flush=True)
    pos = nx.spring_layout(G, k=0.6, iterations=30, seed=42)

    fig, ax = plt.subplots(figsize=(14, 12), facecolor='white')
    
    node_colors = [ccan_color_map[G.nodes[node]['ccan']] for node in G.nodes()]
    degrees = dict(G.degree())
    node_sizes = [(degrees[node] * 150) + 200 for node in G.nodes()]
    node_linewidths = [1.0 + (degrees[node] * 0.3) for node in G.nodes()]

    print("Drawing graph...", flush=True)
    nx.draw_networkx_edges(G, pos, alpha=0.5, width=1.5, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, 
                           edgecolors='black', linewidths=node_linewidths, alpha=0.9, ax=ax)

    labels = {}
    for n in G.nodes():
        if degrees[n] > 2: 
            chrom, coords = n.split(':')
            start = coords.split('-')[0]
            labels[n] = f"{chrom}\n{start}"
            
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight="bold", ax=ax)

    plotted_ccans = set(nx.get_node_attributes(G, 'ccan').values())
    handles = []
    for ccan in unique_ccans:
        if ccan in plotted_ccans:
            count = sum(1 for n in G.nodes() if G.nodes[n]['ccan'] == ccan)
            handles.append(plt.Line2D([0], [0], marker='o', color='w', 
                                      markerfacecolor=ccan_color_map[ccan], 
                                      markersize=12, label=f"CCAN {ccan} ({count} peaks)"))
    
    ax.legend(handles=handles, title="Co-Accessible Networks", loc='upper right', 
              bbox_to_anchor=(1.15, 1), fontsize=10, title_fontsize=12, frameon=False)

    plt.title(f"Chromatin Co-Accessibility Network: {args.prefix}", fontsize=16, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()

    out_path = os.path.join("figures", f"{args.prefix}_network.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Success! Network plot saved to: {out_path}", flush=True)

if __name__ == "__main__":
    main()
