#!/usr/bin/env python3
"""
Co-occurrence analysis: CodY protein alleles vs. Prgg144 binding site marker sequences.

For each pair (CodY allele, Prgg144 sequence), counts the number of strains
that carry both, then outputs a table (CSV) and a heatmap figure.

Usage:
    python cooccurrence_cody_prgg144.py \\
        --cody  results/CodY/allele_info/CodY_Sp_D39/consensus_reference/allele_info_with_strain_name.csv \\
        --prgg  /path/to/marker_region_78_101.csv \\
        --outdir ./cooccurrence_out
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_strain_set(cell, sep=";"):
    """Return a frozenset of stripped strain IDs from a delimited cell."""
    if pd.isna(cell) or str(cell).strip() == "":
        return frozenset()
    return frozenset(x.strip() for x in str(cell).split(sep) if x.strip())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="CodY allele × Prgg144 binding-site co-occurrence analysis"
    )
    parser.add_argument("--cody", required=True,
                        help="allele_info_with_strain_name.csv (CodY)")
    parser.add_argument("--prgg", required=True,
                        help="marker_region CSV with columns: sequence, n_species, genomes")
    parser.add_argument("--outdir", default=".", help="Output directory (created if absent)")
    parser.add_argument("--top_alleles", type=int, default=None,
                        help="Restrict to the top N CodY alleles by frequency (default: all)")
    parser.add_argument("--min_cooccurrence", type=int, default=1,
                        help="Hide allele rows whose max co-occurrence count < this (default: 1)")
    parser.add_argument("--figure_format", default="pdf",
                        choices=["pdf", "png", "svg", "eps"],
                        help="Figure format (default: pdf)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    cody_df = pd.read_csv(args.cody)
    prgg_df = pd.read_csv(args.prgg)

    # Optionally restrict to top N CodY alleles (already sorted by frequency)
    if args.top_alleles:
        cody_df = cody_df.head(args.top_alleles).reset_index(drop=True)

    cody_df["strain_set"] = cody_df["strain_names"].apply(parse_strain_set)
    prgg_df["genome_set"] = prgg_df["genomes"].apply(parse_strain_set)

    # Labels for Prgg144 columns
    prgg_df["label"] = prgg_df["sequence"] + "\n(n=" + prgg_df["n_species"].astype(str) + ")"

    # ------------------------------------------------------------------
    # Build co-occurrence count matrix  (CodY alleles × Prgg sequences)
    # ------------------------------------------------------------------
    n_cody = len(cody_df)
    n_prgg = len(prgg_df)

    counts = np.zeros((n_cody, n_prgg), dtype=int)
    for i, cs in enumerate(cody_df["strain_set"]):
        for j, ps in enumerate(prgg_df["genome_set"]):
            counts[i, j] = len(cs & ps)

    # Filter out allele rows that never co-occur with any Prgg sequence
    row_mask = counts.max(axis=1) >= args.min_cooccurrence
    counts = counts[row_mask]
    cody_df = cody_df[row_mask].reset_index(drop=True)

    # ------------------------------------------------------------------
    # Table output
    # ------------------------------------------------------------------
    table = pd.DataFrame(
        counts,
        index=cody_df["allele_name"],
        columns=prgg_df["sequence"],
    )

    # Add % of CodY allele strains that carry each Prgg sequence (row-normalised)
    for seq in prgg_df["sequence"]:
        table[seq + "_pct_of_cody"] = (
            table[seq] / cody_df["species_count"].values * 100
        ).round(2)

    # Add % of Prgg sequence strains that carry each CodY allele (col-normalised)
    for j, row in prgg_df.iterrows():
        seq = row["sequence"]
        table[seq + "_pct_of_prgg"] = (
            counts[:, j] / row["n_species"] * 100
        ).round(2)

    # Prepend CodY allele metadata columns
    table.insert(0, "cody_species_count", cody_df["species_count"].values)
    table.insert(1, "cody_species_pct", cody_df["species_percentage"].values)
    table.insert(2, "cody_variant_positions",
                 cody_df["variant_positions(position:reference\u2192changed)"].values)

    table_path = os.path.join(args.outdir, "cooccurrence_table.csv")
    table.to_csv(table_path)
    print(f"Table saved: {table_path}")

    # ------------------------------------------------------------------
    # Figure: heatmap (colour = fraction of CodY allele; number = raw count)
    # ------------------------------------------------------------------
    n_rows, n_cols = counts.shape

    # Row-normalise: fraction of each CodY allele's strains carrying each Prgg seq
    row_totals = cody_df["species_count"].values[:, None].astype(float)
    frac = counts / row_totals  # shape (n_rows, n_cols)

    fig_h = max(5, n_rows * 0.32 + 2.5)
    fig_w = max(6, n_cols * 3.0 + 4.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(frac, cmap="Blues", aspect="auto", vmin=0, vmax=1)

    # Annotate each cell with raw count and row-fraction
    for i in range(n_rows):
        for j in range(n_cols):
            raw = counts[i, j]
            pct = frac[i, j] * 100
            text = f"{raw}\n({pct:.0f}%)"
            color = "white" if frac[i, j] > 0.6 else "black"
            ax.text(j, i, text, ha="center", va="center",
                    fontsize=7, color=color, linespacing=1.3)

    # Axes labels
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(prgg_df["label"].tolist(), fontsize=8, ha="center")
    ax.tick_params(axis="x", length=0)

    y_labels = [
        f"{row['allele_name']}  (n={row['species_count']}, {row['species_percentage']:.1f}%)"
        for _, row in cody_df.iterrows()
    ]
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(y_labels, fontsize=7)
    ax.tick_params(axis="y", length=0)

    ax.set_xlabel("Prgg144 CodY-binding site marker sequence", fontsize=10, labelpad=8)
    ax.set_ylabel("CodY protein allele", fontsize=10, labelpad=8)
    ax.set_title(
        "CodY allele × Prgg144 binding-site co-occurrence\n"
        "Color: fraction of CodY allele strains  |  Cell: count (row %)",
        fontsize=10, pad=10,
    )

    # Grid lines between cells
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", length=0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Fraction of CodY allele strains", fontsize=8)
    cbar.ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))

    plt.tight_layout()

    fig_path = os.path.join(args.outdir, f"cooccurrence_figure.{args.figure_format}")
    fig.savefig(fig_path, bbox_inches="tight", dpi=150)
    print(f"Figure saved: {fig_path}")

    # Always also save PNG for quick viewing
    if args.figure_format != "png":
        png_path = os.path.join(args.outdir, "cooccurrence_figure.png")
        fig.savefig(png_path, bbox_inches="tight", dpi=150)
        print(f"Figure saved: {png_path}")

    # ------------------------------------------------------------------
    # Quick summary to stdout
    # ------------------------------------------------------------------
    print("\n--- Co-occurrence count summary ---")
    summary = pd.DataFrame(
        counts,
        index=cody_df["allele_name"],
        columns=prgg_df["sequence"],
    )
    summary.insert(0, "cody_n", cody_df["species_count"].values)
    print(summary.to_string())


if __name__ == "__main__":
    main()
