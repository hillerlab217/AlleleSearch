#!/usr/bin/env python3
"""
Convert allele_info_with_strain_name.csv to a genome × allele binary
presence/absence matrix suitable for upset plots.

Usage:
    python allele_info_to_upset.py allele_info_with_strain_name.csv
    python allele_info_to_upset.py allele_info_with_strain_name.csv -o upset_plot.csv
    python allele_info_to_upset.py allele_info_with_strain_name.csv --min_species_count 2
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def allele_info_to_upset(input_path, output_path, min_species_count=1):
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        sys.exit(f"Error: no rows found in {input_path}")

    # Sort allele rows by their trailing number so columns are ordered 1, 2, 3...
    def allele_number(row):
        m = re.search(r"_(\d+)$", row["allele_name"].strip())
        return int(m.group(1)) if m else float("inf")

    rows.sort(key=allele_number)

    # Apply min_species_count filter
    filtered = []
    for row in rows:
        try:
            count = int(row["species_count"])
        except (KeyError, ValueError):
            print(f"Warning: invalid species_count in row '{row.get('allele_name')}', skipping.")
            continue
        if count >= min_species_count:
            filtered.append(row)
        else:
            print(f"Skipping '{row['allele_name']}' (species_count={count} < {min_species_count})")

    if not filtered:
        sys.exit("Error: no alleles passed the species_count filter.")

    all_allele_names = [row["allele_name"].strip() for row in filtered]

    # Build genome → set of allele names
    genome_to_alleles = defaultdict(set)
    for row in filtered:
        allele_name = row["allele_name"].strip()
        strain_names_raw = row.get("strain_names") or row.get("strain_name", "")
        for genome in strain_names_raw.split(";"):
            genome = genome.strip()
            if genome:
                genome_to_alleles[genome].add(allele_name)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["genome_name"] + all_allele_names)
        writer.writeheader()
        for genome in sorted(genome_to_alleles.keys()):
            row = {"genome_name": genome}
            for allele_name in all_allele_names:
                row[allele_name] = 1 if allele_name in genome_to_alleles[genome] else 0
            writer.writerow(row)

    print(f"Wrote {len(genome_to_alleles)} genomes × {len(all_allele_names)} alleles → {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert allele_info_with_strain_name.csv to an upset plot CSV."
    )
    parser.add_argument("input", help="Path to allele_info_with_strain_name.csv")
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output CSV path (default: upset_plot.csv next to the input file)",
    )
    parser.add_argument(
        "--min_species_count",
        type=int,
        default=1,
        help="Minimum species count for an allele to be included (default: 1)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"Error: input file not found: {input_path}")

    output_path = Path(args.output) if args.output else input_path.parent / "upset_plot.csv"

    allele_info_to_upset(input_path, output_path, args.min_species_count)


if __name__ == "__main__":
    main()
