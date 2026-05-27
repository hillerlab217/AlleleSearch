#!/usr/bin/env python3
"""
Add per-position nucleotide columns for CodY_Binding_Site alleles into a per-strain CSV.

For each strain, looks up their CodY_Binding_Site allele number, then reads the
nucleotide at the specified positions (1-based) from the allele sequence.
Strains with CodY_Binding_Site == 0 get 0 in all position columns.

Usage:
    python add_binding_site_positions.py \
        --per_strain_sheet data/allele_each_strain_sp.csv \
        --allele_info /path/to/allele_info_no_strain_name.csv \
        --source_column CodY_Binding_Site \
        --positions 3 5 13 24
"""

import argparse
import csv
import re
import sys


def load_allele_sequences(allele_info_path, source_column_prefix):
    """
    Parse allele_info_no_strain_name.csv and return a dict:
        allele_number (int) -> sequence (str)
    """
    allele_seq = {}
    with open(allele_info_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            allele_name = row["allele_name"].strip()
            m = re.match(r"^(.+)_(\d+)$", allele_name)
            if not m:
                print(f"Warning: cannot parse allele number from '{allele_name}', skipping.")
                continue
            number = int(m.group(2))
            allele_seq[number] = row["allele_sequence"].strip()
    return allele_seq


def main():
    parser = argparse.ArgumentParser(
        description="Add per-position nucleotide columns for a binding site allele column."
    )
    parser.add_argument(
        "--per_strain_sheet",
        required=True,
        help="Path to the per-strain CSV (modified in-place).",
    )
    parser.add_argument(
        "--allele_info",
        required=True,
        help="Path to allele_info_no_strain_name.csv with allele sequences.",
    )
    parser.add_argument(
        "--source_column",
        required=True,
        help="Column in per-strain sheet holding the allele number (e.g. CodY_Binding_Site).",
    )
    parser.add_argument(
        "--positions",
        type=int,
        nargs="+",
        required=True,
        help="1-based positions in the allele sequence to extract (e.g. 3 5 13 24).",
    )
    args = parser.parse_args()

    allele_seq = load_allele_sequences(args.allele_info, args.source_column)

    # Build new column names: {source_column}_P{pos}
    new_cols = [f"{args.source_column}_P{p}" for p in args.positions]

    # Validate all sequences cover the requested positions
    max_pos = max(args.positions)
    for num, seq in allele_seq.items():
        if len(seq) < max_pos:
            print(
                f"Warning: allele {num} sequence length {len(seq)} < max position {max_pos}."
            )

    # Read existing per-strain sheet
    with open(args.per_strain_sheet, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    if args.source_column not in fieldnames:
        sys.exit(f"Error: column '{args.source_column}' not found in {args.per_strain_sheet}")

    # Add new columns (overwrite if already present)
    for col in new_cols:
        if col not in fieldnames:
            fieldnames.append(col)

    # Fill in values
    for row in rows:
        allele_num_raw = row[args.source_column]
        try:
            allele_num = int(allele_num_raw)
        except (ValueError, TypeError):
            allele_num = 0

        if allele_num == 0:
            for col in new_cols:
                row[col] = 0
            continue

        seq = allele_seq.get(allele_num)
        if seq is None:
            print(f"Warning: allele {allele_num} not found in allele_info, setting 0.")
            for col in new_cols:
                row[col] = 0
            continue

        for pos, col in zip(args.positions, new_cols):
            row[col] = seq[pos - 1]  # convert 1-based to 0-based

    # Write back in-place
    with open(args.per_strain_sheet, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Added columns {new_cols} to {args.per_strain_sheet}")


if __name__ == "__main__":
    main()
