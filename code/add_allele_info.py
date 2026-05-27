#!/usr/bin/env python3
"""
Add allele information from an allele_info_with_strain_name.csv into a per-strain CSV.

Usage:
    python add_allele_info.py \
        --per_strain_sheet allele_search/results/allele_each_strain_sp.csv \
        allele_search/results/CodY/allele_info/CodY_Sp_D39/self_reference/allele_info_with_strain_name.csv \
        [more allele_info files ...]

The new column header is derived from the allele_name prefix (e.g. CodY_Sp_D39_1 -> CodY_Sp_D39).
Strains in allele N get value N; strains absent from all alleles get 0.
The modified CSV is written in-place.
"""

import argparse
import csv
import re
import sys


def screen_min_species_count(rows, min_species_count):
    """
    Filter allele rows by minimum species count.

    Args:
        rows (list[dict]): rows from allele_info CSV (each must have 'species_count')
        min_species_count (int): minimum species count; rows below this are excluded

    Returns:
        list[dict]: rows whose species_count >= min_species_count
    """
    kept = []
    for row in rows:
        try:
            count = int(row["species_count"])
        except (KeyError, ValueError):
            print(f"Warning: missing or invalid 'species_count' in row {row}, skipping.")
            continue
        if count >= min_species_count:
            kept.append(row)
        else:
            print(
                f"Filtered out allele '{row.get('allele_name', '?')}' "
                f"(species_count={count} < {min_species_count})"
            )
    return kept


def parse_allele_info(allele_info_path, min_species_count=1):
    """
    Returns:
        col_header (str): e.g. 'CodY_Sp_D39'
        strain_to_allele (dict): strain_name -> allele number (int)
    """
    strain_to_allele = {}
    col_header = None

    with open(allele_info_path, newline="") as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)

    filtered_rows = screen_min_species_count(all_rows, min_species_count)

    for row in filtered_rows:
        allele_name = row["allele_name"].strip()

        # Derive header and allele number: last '_N' segment is the number
        m = re.match(r"^(.+)_(\d+)$", allele_name)
        if not m:
            print(f"Warning: cannot parse allele number from '{allele_name}', skipping.")
            continue

        prefix, number = m.group(1), int(m.group(2))

        if col_header is None:
            col_header = prefix
        elif col_header != prefix:
            print(
                f"Warning: mixed prefixes in {allele_info_path} "
                f"('{col_header}' vs '{prefix}'). Using first."
            )

        strain_names_raw = row.get("strain_names") or row.get("strain_name", "")
        for strain in strain_names_raw.split(";"):
            strain = strain.strip()
            if strain:
                strain_to_allele[strain] = number

    if col_header is None:
        sys.exit(f"Error: no valid allele rows found in {allele_info_path}")

    return col_header, strain_to_allele


def main():
    parser = argparse.ArgumentParser(description="Add allele columns to a per-strain CSV.")
    parser.add_argument(
        "--per_strain_sheet",
        required=True,
        help="Path to the per-strain CSV (modified in-place).",
    )
    parser.add_argument(
        "allele_info_files",
        nargs="+",
        help="One or more allele_info_with_strain_name.csv files.",
    )
    parser.add_argument(
        "--column_name",
        default=None,
        help="Custom name for the new column. Only valid when providing a single allele info file.",
    )
    parser.add_argument(
        "--min_species_count",
        type=int,
        default=1,
        help="Minimum species count for an allele to be included (default: 1).",
    )
    args = parser.parse_args()

    if args.column_name and len(args.allele_info_files) > 1:
        sys.exit("Error: --column_name can only be used with a single allele info file.")

    # Read existing per-strain sheet
    with open(args.per_strain_sheet, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    # Process each allele info file
    for allele_info_path in args.allele_info_files:
        col_header, strain_to_allele = parse_allele_info(allele_info_path, args.min_species_count)
        if args.column_name:
            col_header = args.column_name

        if col_header in fieldnames:
            print(f"Column '{col_header}' already exists — overwriting.")
        else:
            fieldnames.append(col_header)

        for row in rows:
            strain = row["strain_name"].strip()
            row[col_header] = strain_to_allele.get(strain, 0)

        print(f"Added column '{col_header}' from {allele_info_path}")

    # Write back in-place
    with open(args.per_strain_sheet, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {args.per_strain_sheet}")


if __name__ == "__main__":
    main()
