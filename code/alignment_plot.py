#!/usr/bin/env python3
"""
extract_variant_aas_excel.py

Build an Excel alignment table.

Options
-------
• With the extra argument  gap  → keep every position from min to max,
  leaving non-variant columns blank.
• Without that argument       → include only variant-site columns.

Usage
-----
python extract_variant_aas_excel.py INPUT.csv OUTPUT.xlsx MIN_PERCENTAGE [gap]
"""

import sys
import re
import pandas as pd


# -------------------------- helper ---------------------------------------
def extract_positions(variant_str: str) -> list[int]:
    """'166:D→E; 469:V→I' → [166, 469]"""
    if pd.isna(variant_str) or not variant_str.strip():
        return []
    return [
        int(m.group(1))
        for part in re.split(r"[;,]", variant_str)
        if (m := re.match(r"\s*(\d+)", part))
    ]


# --------------------------- main ----------------------------------------
def main(input_csv: str, output_xlsx: str, min_pct: float, keep_gaps: bool) -> None:
    df = pd.read_csv(input_csv)

    PCT_COL, VARIANT_COL, SEQ_COL, NAME_COL = (
        "species_percentage",
        "variant_positions(position:reference→changed)",
        "allele_sequence",
        "allele_name",
    )

    # ---- checks ---------------------------------------------------------
    for col in (PCT_COL, VARIANT_COL, SEQ_COL, NAME_COL):
        if col not in df.columns:
            sys.exit(f"Column '{col}' not found in {input_csv}")

    # ---- filter rows ----------------------------------------------------
    df = df[df[PCT_COL] >= min_pct].copy()
    if df.empty:
        print(f"No rows meet ≥{min_pct}% — nothing written.")
        return

    # ---- collect variant positions -------------------------------------
    variant_positions = set()
    for s in df[VARIANT_COL].fillna(""):
        variant_positions.update(extract_positions(s))

    if not variant_positions:
        print("No variant positions found — nothing written.")
        return

    # ---- decide column set ---------------------------------------------
    if keep_gaps:
        columns_to_keep = list(range(min(variant_positions), max(variant_positions) + 1))
    else:
        columns_to_keep = sorted(variant_positions)

    # ---- build table ----------------------------------------------------
    rows_out = []
    for _, row in df.iterrows():
        seq = row[SEQ_COL]
        aa_cells = [
            seq[p - 1] if (p in variant_positions and p - 1 < len(seq)) else ""
            for p in columns_to_keep
        ]
        rows_out.append([row[NAME_COL], *aa_cells])

    headers = ["allele_name"] + [
        str(p) if p in variant_positions else "" for p in columns_to_keep
    ]

    pd.DataFrame(rows_out, columns=headers).to_excel(
        output_xlsx, index=False, sheet_name="alignment"
    )
    print(f"✓ Wrote {output_xlsx}")


# --------------------------- CLI -----------------------------------------
if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        sys.exit(
            "Usage: python extract_variant_aas_excel.py INPUT.csv OUTPUT.xlsx "
            "MIN_PERCENTAGE [gap]"
        )

    in_csv, out_xlsx, pct_str = sys.argv[1:4]
    try:
        pct_threshold = float(pct_str)
    except ValueError:
        sys.exit(f"MIN_PERCENTAGE must be numeric, got '{pct_str}'")

    keep_gaps = (
        len(sys.argv) == 5
        and sys.argv[4].lower() in {"gap", "gaps", "true", "yes", "1"}
    )

    main(in_csv, out_xlsx, pct_threshold, keep_gaps)
