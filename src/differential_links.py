#!/usr/bin/env python3

import sys
import argparse
import os
import pandas as pd
import itertools

# Force output to stderr so it NEVER buffers or hides
print("========================================", file=sys.stderr, flush=True)
print("SCRIPT STARTED. READING ARGUMENTS...", file=sys.stderr, flush=True)
print("========================================", file=sys.stderr, flush=True)


def peak_sort_key(peak):
    chrom, coords = peak.split(':')
    start, end = coords.split('-')
    return chrom, int(start), int(end)


def main():
    parser = argparse.ArgumentParser(description="Calculate ALL pairwise differential co-accessibility links with summary.")
    parser.add_argument("--csv", nargs='+', required=True, help="List of links CSV files.")
    parser.add_argument("--prefix", required=True, help="Prefix for the output files.")
    parser.add_argument("--cutoff", type=float, required=True, help="Existence threshold on CIRCE score. Links with score <= cutoff are considered absent in that condition.")
    args = parser.parse_args()

    print(f"Arguments parsed successfully: {args}", file=sys.stderr, flush=True)

    cwd = os.getcwd()
    print(f"Working Directory: {cwd}", file=sys.stderr, flush=True)

    dfs = {}

    # 1. Load all files
    for f in args.csv:
        cond_name = os.path.basename(f).replace("circe_", "").replace("_coaccessibility_links.csv", "").replace(".csv", "")
        print(f"Loading {cond_name} from {f}...", file=sys.stderr, flush=True)

        if not os.path.exists(f):
            print(f"ERROR: FILE DOES NOT EXIST: {f}", file=sys.stderr, flush=True)
            sys.exit(1)

        df = pd.read_csv(f)
        print(f"  -> Loaded {len(df):,} rows from {f}", file=sys.stderr, flush=True)

        score_col = 'score' if 'score' in df.columns else 'coaccess'
        df = df[['Peak1', 'Peak2', score_col]].rename(columns={score_col: cond_name})

        # Sort peaks using numeric genomic coordinates
        ordered_peaks = df.apply(
            lambda row: sorted((row['Peak1'], row['Peak2']), key=peak_sort_key),
            axis=1
        )
        df['Peak1'] = [p[0] for p in ordered_peaks]
        df['Peak2'] = [p[1] for p in ordered_peaks]

        # Report duplicate links before removing them
        duplicate_count = df.duplicated(subset=['Peak1', 'Peak2']).sum()

        if duplicate_count > 0:
            print(
                f"  -> WARNING: Found {duplicate_count:,} duplicate links in {cond_name}. "
                f"Keeping the first occurrence of each link.",
                file=sys.stderr,
                flush=True
            )

        df = df.drop_duplicates(subset=['Peak1', 'Peak2'])

        # Apply existence cutoff: links at or below cutoff are considered absent
        before_cutoff = len(df)
        df = df[df[cond_name] > args.cutoff]
        print(
            f"  -> Existence cutoff {args.cutoff}: kept {len(df):,} of {before_cutoff:,} links "
            f"({before_cutoff - len(df):,} removed as below cutoff).",
            file=sys.stderr,
            flush=True
        )

        dfs[cond_name] = df

    cond_names = list(dfs.keys())
    print(f"Files loaded successfully. Conditions: {cond_names}", file=sys.stderr, flush=True)

    # 2. Merge all raw data
    print("Merging dataframes (this may take a minute)...", file=sys.stderr, flush=True)

    merged = dfs[cond_names[0]]

    for name in cond_names[1:]:
        print(f"  Merging {name}...", file=sys.stderr, flush=True)
        merged = pd.merge(merged, dfs[name], on=['Peak1', 'Peak2'], how='outer')

    print(f"Merge complete. Total unique links: {len(merged):,}", file=sys.stderr, flush=True)

    # 3. Calculate Status for EVERY possible pair (based on existence / null-ness)
    pairs = list(itertools.combinations(cond_names, 2))

    for (A, B) in pairs:
        comp_name = f"{A}_vs_{B}"
        status_col = f"Status_{comp_name}"

        present_A = merged[A].notna()
        present_B = merged[B].notna()

        def assign_status(a_present, b_present, b_name):
            if a_present and b_present:
                return "Shared"
            elif a_present and not b_present:
                return f"Lost_in_{b_name}"
            elif not a_present and b_present:
                return f"Gained_in_{b_name}"

        merged[status_col] = [
            assign_status(a, b, B)
            for a, b in zip(present_A, present_B)
        ]

    # --- AGGREGATE SUMMARY PRINT BLOCK ---
    print("\n" + "=" * 50, file=sys.stderr, flush=True)
    print("PAIRWISE DIFFERENTIAL SUMMARY STATISTICS:", file=sys.stderr, flush=True)
    print("=" * 50, file=sys.stderr, flush=True)

    for (A, B) in pairs:
        comp_name = f"{A}_vs_{B}"
        status_col = f"Status_{comp_name}"
        counts = merged[status_col].value_counts()

        # Fixed print order: Shared, Lost_in_B, Gained_in_B
        ordered_keys = [
            "Shared",
            f"Lost_in_{B}",
            f"Gained_in_{B}",
        ]

        print(f"\nComparison: {A} vs {B}", file=sys.stderr, flush=True)
        for status in ordered_keys:
            count = counts.get(status, 0)
            print(f"  - {status}: {count:,}", file=sys.stderr, flush=True)

    print("\n" + "=" * 50, file=sys.stderr, flush=True)
    # -------------------------------------

    # 4. Save CSV in CURRENT DIRECTORY
    out_csv = os.path.join(cwd, f"{args.prefix}_all_pairwise_differential_links.csv")

    print(f"Saving CSV to: {out_csv}", file=sys.stderr, flush=True)

    merged.to_csv(out_csv, index=False)

    print(f"SUCCESS! Differential statistics table saved to {out_csv}.", file=sys.stderr, flush=True)


if __name__ == "__main__":
    try:
        main()

        print("========================================", file=sys.stderr, flush=True)
        print("SCRIPT COMPLETED SUCCESSFULLY.", file=sys.stderr, flush=True)
        print("========================================", file=sys.stderr, flush=True)

    except Exception as e:
        import traceback

        print("\n========================================", file=sys.stderr, flush=True)
        print("FATAL ERROR ENCOUNTERED:", file=sys.stderr, flush=True)
        traceback.print_exc(file=sys.stderr)
        print("========================================", file=sys.stderr, flush=True)

        sys.exit(1)

