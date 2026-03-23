#!/usr/bin/env python3
"""Filter a BLAST CSV by E-value, sequence identity, and alignment length."""
import csv
import argparse
from collections import defaultdict
from Bio import SeqIO


def compute_min_length(query_fasta, length_percentage):
    """Return minimum alignment length: total query nt * length_percentage / 100."""
    total = sum(len(rec.seq) for rec in SeqIO.parse(query_fasta, 'fasta'))
    return int(total * length_percentage / 100)


def filter_blast(input_file, output_file, max_evalue, sequence_identity, min_length):
    strain_hits = defaultdict(int)

    with open(input_file) as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)
        header = next(reader)
        writer.writerow(header)
        ei = header.index('E-value')
        li = header.index('Alignment Length')
        ii = header.index('ID')
        pi = (header.index('Sequence Identity')
              if 'Sequence Identity' in header
              else header.index('Alignment Percentage'))
        for row in reader:
            if (float(row[ei]) < max_evalue
                    and float(row[pi]) >= sequence_identity
                    and int(row[li]) >= min_length):
                writer.writerow(row)
                strain = row[ii].split('_vs_', 1)[-1]
                strain_hits[strain] += 1

    duplicates = [(s, n) for s, n in strain_hits.items() if n > 1]
    if duplicates:
        print(f"WARNING: {len(duplicates)} strain(s) have multiple hits passing screening "
              f"(possible gene duplication):")
        for strain, count in sorted(duplicates):
            print(f"  {strain}: {count} hits")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Filter BLAST CSV by E-value, sequence identity, and alignment length.')
    parser.add_argument('--input', required=True, help='Input BLAST CSV file')
    parser.add_argument('--output', required=True, help='Output filtered CSV file')
    parser.add_argument('--max_evalue', type=float, required=True,
                        help='Maximum E-value threshold (exclusive)')
    parser.add_argument('--sequence_identity', type=float, required=True,
                        help='Minimum sequence identity %%')
    parser.add_argument('--query_fasta', required=True,
                        help='Query FASTA file (used to compute minimum alignment length)')
    parser.add_argument('--length_percentage', type=float, required=True,
                        help='Minimum alignment length as %% of query length')
    args = parser.parse_args()

    min_length = compute_min_length(args.query_fasta, args.length_percentage)
    filter_blast(args.input, args.output, args.max_evalue,
                 args.sequence_identity, min_length)
