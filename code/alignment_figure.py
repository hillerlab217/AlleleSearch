#!/usr/bin/env python3
"""
alignment_figure.py

Draw a coloured amino-acid variant alignment figure directly from
the allele-info CSV (no Excel intermediate).

Amino acids are coloured by physicochemical group (same palette as
the original Excel workflow).  Only variant-position columns are
shown unless ``--gap`` is given, which fills in every position from
the first to the last variant site.

Usage
-----
python alignment_figure.py INPUT.csv OUTPUT.{pdf,png,svg} MIN_PERCENTAGE [--gap]

Examples
--------
python alignment_figure.py allele_info.csv alignment.pdf 1.0
python alignment_figure.py allele_info.csv alignment.png 1.0 --gap
"""

import sys
import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgb


# ── amino-acid colour palette (physicochemical groups) ──────────────
AA_HEX = {
    # Hydrophobic / aliphatic
    "A": "#f4b183", "G": "#f4b183", "I": "#f4b183", "L": "#f4b183",
    "M": "#f4b183", "V": "#f4b183", "P": "#f4b183",
    # Positive charge
    "R": "#4f81bd", "K": "#4f81bd", "H": "#4f81bd",
    # Polar / uncharged
    "N": "#9bbb59", "C": "#9bbb59", "Q": "#9bbb59", "S": "#9bbb59", "T": "#9bbb59",
    # Negative charge
    "D": "#e26b6b", "E": "#e26b6b",
    # Aromatic
    "F": "#fff475", "W": "#fff475", "Y": "#fff475",
}

LEGEND_GROUPS = {
    "Hydrophobic (A G I L M V P)": "#f4b183",
    "Positive (R K H)":           "#4f81bd",
    "Polar (N C Q S T)":          "#9bbb59",
    "Negative (D E)":             "#e26b6b",
    "Aromatic (F W Y)":           "#fff475",
}

GAP_RGB   = to_rgb("#eeeeee")
BLANK_RGB = to_rgb("#ffffff")


# ── helpers ─────────────────────────────────────────────────────────
def extract_positions(variant_str: str) -> list[int]:
    """'166:D→E; 469:V→I' → [166, 469]"""
    if pd.isna(variant_str) or not variant_str.strip():
        return []
    return [
        int(m.group(1))
        for part in re.split(r"[;,]", variant_str)
        if (m := re.match(r"\s*(\d+)", part))
    ]


# ── data extraction ─────────────────────────────────────────────────
def build_alignment_data(input_csv: str, min_pct: float, keep_gaps: bool):
    """Return (row_labels, col_positions, variant_positions, grid)."""
    df = pd.read_csv(input_csv)

    PCT_COL     = "species_percentage"
    VARIANT_COL = "variant_positions(position:reference→changed)"
    SEQ_COL     = "allele_sequence"
    NAME_COL    = "allele_name"

    for col in (PCT_COL, VARIANT_COL, SEQ_COL, NAME_COL):
        if col not in df.columns:
            sys.exit(f"Column '{col}' not found in {input_csv}")

    df = df[df[PCT_COL] >= min_pct].copy()
    if df.empty:
        sys.exit(f"No rows meet ≥{min_pct}% — nothing to draw.")

    variant_positions = set()
    for s in df[VARIANT_COL].fillna(""):
        variant_positions.update(extract_positions(s))
    if not variant_positions:
        sys.exit("No variant positions found — nothing to draw.")

    if keep_gaps:
        columns = list(range(min(variant_positions), max(variant_positions) + 1))
    else:
        columns = sorted(variant_positions)

    row_labels = []
    grid = []
    for _, row in df.iterrows():
        seq = row[SEQ_COL]
        cells = [
            seq[p - 1] if (p in variant_positions and p - 1 < len(seq)) else ""
            for p in columns
        ]
        row_labels.append(str(row[NAME_COL]))
        grid.append(cells)

    return row_labels, columns, variant_positions, grid


# ── drawing ─────────────────────────────────────────────────────────
def draw_alignment(row_labels, columns, variant_positions, grid, output_path):
    """Render the alignment grid using imshow for pixel-perfect lines."""
    n_rows = len(row_labels)
    n_cols = len(columns)
    is_variant = [c in variant_positions for c in columns]

    # Build RGB colour matrix  (n_rows × n_cols × 3)
    img = np.ones((n_rows, n_cols, 3))
    for i in range(n_rows):
        for j in range(n_cols):
            if not is_variant[j]:
                img[i, j] = GAP_RGB
            else:
                aa = grid[i][j].upper() if grid[i][j] else ""
                img[i, j] = to_rgb(AA_HEX[aa]) if aa in AA_HEX else BLANK_RGB

    # ── figure sizing ───────────────────────────────────────────────
    # In gap mode, make non-variant columns narrower via column width ratios
    if any(not v for v in is_variant):
        # gap mode: variant cols = 1.0, gap cols = 0.2 relative width
        col_w = [1.0 if v else 0.2 for v in is_variant]
    else:
        col_w = [1.0] * n_cols

    total_col_units = sum(col_w)

    # Scale so variant cells are ~0.3 inches wide
    cell_inch = 0.30
    grid_w_inch = cell_inch * total_col_units
    # cap max width at 40 inches for sanity
    if grid_w_inch > 40:
        cell_inch = 40 / total_col_units
        grid_w_inch = 40

    cell_h_inch = 0.40
    label_margin = 1.6   # inches for allele names on left
    header_margin = 1.0  # inches for rotated position numbers on top
    legend_h = 0.5

    fig_w = label_margin + grid_w_inch + 0.3
    fig_h = legend_h + header_margin + cell_h_inch * n_rows + 0.2

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    # axes for the grid image
    ax_left   = label_margin / fig_w
    ax_bottom = 0.2 / fig_h
    ax_width  = grid_w_inch / fig_w
    ax_height = (cell_h_inch * n_rows) / fig_h

    ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])

    # ── build x-edge array for pcolormesh (variable-width columns) ──
    x_edges = [0.0]
    for w in col_w:
        x_edges.append(x_edges[-1] + w)
    x_edges = np.array(x_edges)
    y_edges = np.arange(n_rows + 1, dtype=float)

    # pcolormesh needs (M+1, N+1) edges and (M, N) shaped color data
    # Flip rows so row 0 is at top
    img_flipped = img[::-1]

    # Reshape colour data for pcolormesh: flatten RGB into facecolors
    # Use pcolormesh with shading='flat'
    ax.pcolormesh(
        x_edges, y_edges, img_flipped,
        shading="flat",
        edgecolors="#dddddd",
        linewidth=0.3,
    )

    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    ax.invert_yaxis()

    # ── AA text labels inside cells ─────────────────────────────────
    font_aa = 7 if n_cols > 80 else 8
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = np.arange(n_rows) + 0.5

    for i in range(n_rows):
        for j in range(n_cols):
            if is_variant[j]:
                aa = grid[i][j].upper() if grid[i][j] else ""
                if aa:
                    ax.text(
                        x_centers[j], y_centers[i], aa,
                        ha="center", va="center",
                        fontsize=font_aa,
                    )

    # ── y-axis: allele names ────────────────────────────────────────
    ax.set_yticks(y_centers)
    ax.set_yticklabels(row_labels, fontsize=font_aa)
    ax.tick_params(axis="y", length=0, pad=4)

    # ── x-axis: position numbers on top ─────────────────────────────
    # Only label variant columns
    variant_ticks = [x_centers[j] for j in range(n_cols) if is_variant[j]]
    variant_labels = [str(columns[j]) for j in range(n_cols) if is_variant[j]]

    ax.xaxis.set_ticks_position("top")
    ax.set_xticks(variant_ticks)
    ax.set_xticklabels(
        variant_labels,
        fontsize=font_aa - 1, fontweight="bold",
        rotation=90, ha="center",
    )
    ax.tick_params(axis="x", length=0, pad=2)

    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    # ── legend above the axes ───────────────────────────────────────
    legend_handles = [
        mpatches.Patch(facecolor=c, edgecolor="grey", linewidth=0.5, label=lab)
        for lab, c in LEGEND_GROUPS.items()
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.01, 0.99),
        ncol=min(len(LEGEND_GROUPS), 3),
        fontsize=6, frameon=False,
        handlelength=1.2, handletextpad=0.4, columnspacing=1.0,
    )

    # ── save ────────────────────────────────────────────────────────
    suffix = Path(output_path).suffix.lower()
    dpi = 300 if suffix == ".png" else 150
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {output_path}")


# ── CLI ─────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Draw a coloured amino-acid variant alignment figure.",
    )
    parser.add_argument("input_csv", help="Allele-info CSV file")
    parser.add_argument("output", help="Output figure path (.pdf, .png, or .svg)")
    parser.add_argument("min_percentage", type=float,
                        help="Minimum species %% to include")
    parser.add_argument("--gap", action="store_true",
                        help="Show all positions between first and last variant")
    args = parser.parse_args()

    row_labels, columns, variant_positions, grid = build_alignment_data(
        args.input_csv, args.min_percentage, args.gap,
    )
    draw_alignment(row_labels, columns, variant_positions, grid, args.output)


if __name__ == "__main__":
    main()
