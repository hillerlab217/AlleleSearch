#!/usr/bin/env python3
"""
Analyze blast summary filtered CSV files to:
1. Count unique strains per region file and list any replicates (strains with multiple hits)
2. Group unique sequences and count how many strains share each unique sequence

Usage:
    python strain_sequence_analysis.py <input_directory> [output_directory]

Example:
    python strain_sequence_analysis.py results/s_p_insertion/blast_summary_filtered_98
"""

import csv
import os
import sys
import hashlib
from collections import defaultdict


def extract_strain_id(row_id):
    """Extract strain ID from the blast result ID (everything after '_vs_')."""
    parts = row_id.split("_vs_")
    if len(parts) >= 2:
        return "_".join(parts[1:])
    return row_id


def analyze_file(filepath):
    """Parse a single CSV file and return strain->sequence mapping."""
    strain_sequences = defaultdict(list)

    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            strain_id = extract_strain_id(row["ID"])
            sequence = row["Sequence"].strip()
            strain_sequences[strain_id].append(sequence)

    return strain_sequences


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input_directory> [output_directory]")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else input_dir

    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    csv_files = sorted(
        [f for f in os.listdir(input_dir) if f.endswith(".csv")]
    )

    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        sys.exit(1)

    # --- Output file 1: Strain summary + replicates ---
    strain_report_path = os.path.join(output_dir, "strain_summary.txt")
    with open(strain_report_path, "w") as report:
        report.write("=" * 70 + "\n")
        report.write("STRAIN SUMMARY REPORT\n")
        report.write("=" * 70 + "\n\n")

        for csv_file in csv_files:
            filepath = os.path.join(input_dir, csv_file)
            strain_sequences = analyze_file(filepath)

            region_name = os.path.splitext(csv_file)[0]
            report.write(f"--- {region_name} ---\n")
            report.write(f"Total hits: {sum(len(v) for v in strain_sequences.values())}\n")
            report.write(f"Unique strains: {len(strain_sequences)}\n")

            # Find replicates (strains with more than one hit)
            replicates = {
                strain: seqs
                for strain, seqs in strain_sequences.items()
                if len(seqs) > 1
            }

            if replicates:
                report.write(f"Replicate strains (multiple hits): {len(replicates)}\n")
                for strain, seqs in sorted(replicates.items()):
                    same_seq = len(set(seqs)) == 1
                    report.write(
                        f"  {strain}: {len(seqs)} hits"
                        f" ({'identical sequences' if same_seq else 'different sequences'})\n"
                    )
            else:
                report.write("Replicate strains: None\n")

            report.write("\n")

    print(f"Strain summary written to: {strain_report_path}")

    # --- Output file 2: Unique sequences and their strain counts ---
    seq_report_path = os.path.join(output_dir, "unique_sequences.csv")
    with open(seq_report_path, "w", newline="") as seq_out:
        writer = csv.writer(seq_out)
        writer.writerow([
            "Region", "Sequence_ID", "Strain_Count", "Strains", "Sequence"
        ])

        for csv_file in csv_files:
            filepath = os.path.join(input_dir, csv_file)
            strain_sequences = analyze_file(filepath)
            region_name = os.path.splitext(csv_file)[0]

            # Group strains by their sequence
            seq_to_strains = defaultdict(list)
            for strain, seqs in strain_sequences.items():
                for seq in set(seqs):  # use set to avoid counting same strain twice for identical seqs
                    seq_to_strains[seq].append(strain)

            # Sort by number of strains (descending)
            for idx, (seq, strains) in enumerate(
                sorted(seq_to_strains.items(), key=lambda x: -len(x[1])), 1
            ):
                seq_hash = hashlib.md5(seq.encode()).hexdigest()[:8]
                seq_id = f"{region_name}_seq_{idx}_{seq_hash}"
                writer.writerow([
                    region_name,
                    seq_id,
                    len(strains),
                    "; ".join(sorted(strains)),
                    seq,
                ])

    print(f"Unique sequences written to: {seq_report_path}")

    # --- Also print a concise summary to stdout ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    for csv_file in csv_files:
        filepath = os.path.join(input_dir, csv_file)
        strain_sequences = analyze_file(filepath)
        region_name = os.path.splitext(csv_file)[0]

        seq_to_strains = defaultdict(list)
        for strain, seqs in strain_sequences.items():
            for seq in set(seqs):
                seq_to_strains[seq].append(strain)

        print(f"\n{region_name}:")
        print(f"  Total hits:        {sum(len(v) for v in strain_sequences.values())}")
        print(f"  Unique strains:    {len(strain_sequences)}")
        print(f"  Unique sequences:  {len(seq_to_strains)}")

        replicates = [s for s, v in strain_sequences.items() if len(v) > 1]
        if replicates:
            print(f"  Replicate strains: {len(replicates)} ({', '.join(sorted(replicates))})")

        # Show top sequences by strain count
        top_seqs = sorted(seq_to_strains.items(), key=lambda x: -len(x[1]))[:5]
        print(f"  Top sequences by strain count:")
        for i, (seq, strains) in enumerate(top_seqs, 1):
            seq_preview = seq[:40] + "..." if len(seq) > 40 else seq
            print(f"    #{i}: {len(strains)} strains | {seq_preview}")


if __name__ == "__main__":
    main()
