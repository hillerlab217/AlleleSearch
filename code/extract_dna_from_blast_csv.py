"""
Extract DNA sequences from genome files using BLAST subject coordinates.

When tblastn is used, the BLAST CSV contains protein sequences but with
nucleotide coordinates. This script replaces the protein sequences with
the actual DNA extracted from the genome at those coordinates.

Usage:
    python extract_dna_from_blast_csv.py <input_csv> <output_csv> <genome_dir>
"""

import csv
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq


def find_genome_file(genome_dir, genome_name):
    """Find the genome FASTA file for a given genome name."""
    for ext in ['.fasta', '.fa', '.fna', '.ffn']:
        path = os.path.join(genome_dir, genome_name + ext)
        if os.path.exists(path):
            return path
    return None


def load_genome(genome_file):
    """Load all contigs from a genome FASTA into a dict keyed by ID."""
    contigs = {}
    for record in SeqIO.parse(genome_file, 'fasta'):
        contigs[record.id] = record
    return contigs


def find_contig(contigs, contig_name):
    """Look up a contig with ordered fallbacks for real-world ID mismatches.

    Fallback order:
    1. Exact match
    2. Case-insensitive exact match
    3. Version-stripped match (e.g. NZ_ABC.1 matches NZ_ABC.2)
    4. Bidirectional substring match (last resort)
    """
    if contig_name in contigs:
        return contigs[contig_name]
    lower = contig_name.lower()
    for cid, rec in contigs.items():
        if cid.lower() == lower:
            return rec
    stripped = contig_name.rsplit('.', 1)[0]
    for cid, rec in contigs.items():
        if cid.rsplit('.', 1)[0] == stripped:
            return rec
    for cid, rec in contigs.items():
        if contig_name in cid or cid in contig_name:
            return rec
    return None


def extract_dna(contigs, contig_name, start, end):
    """Extract DNA from a contig at the given 1-based coordinates.

    If start < end: forward strand, extract seq[start-1:end]
    If start > end: reverse strand, extract seq[end-1:start] then reverse complement
    """
    record = find_contig(contigs, contig_name)

    if record is None:
        return None

    if start <= end:
        dna = str(record.seq[start - 1:end])
    else:
        dna_raw = str(record.seq[end - 1:start])
        dna = str(Seq(dna_raw).reverse_complement())

    return dna


def main():
    parser = argparse.ArgumentParser(
        description="Replace protein sequences in BLAST CSV with DNA from genomes."
    )
    parser.add_argument("input_csv", help="Filtered BLAST CSV with protein sequences")
    parser.add_argument("output_csv", help="Output CSV with DNA sequences")
    parser.add_argument("genome_dir", help="Directory containing genome FASTA files")
    args = parser.parse_args()

    # Cache loaded genomes to avoid re-reading
    genome_cache = {}

    with open(args.input_csv) as infile, open(args.output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()

        for row in reader:
            blast_id = row['ID']
            contig = row['Contig']
            start = int(row['Subject Start'])
            end = int(row['Subject End'])

            # Extract genome name from ID: {query}_vs_{genome_name}
            genome_name = blast_id.split("_vs_", 1)[-1]

            # Load genome if not cached
            if genome_name not in genome_cache:
                genome_file = find_genome_file(args.genome_dir, genome_name)
                if genome_file:
                    genome_cache[genome_name] = load_genome(genome_file)
                else:
                    print(f"WARNING: Genome file not found for {genome_name}, skipping row")
                    continue

            contigs = genome_cache[genome_name]
            dna = extract_dna(contigs, contig, start, end)

            if dna is None:
                print(f"WARNING: Contig {contig} not found in {genome_name}, skipping row")
                continue

            # Replace protein sequence with DNA and update alignment length
            row['Sequence'] = dna
            row['Alignment Length'] = str(len(dna))
            writer.writerow(row)

    print(f"DNA extraction complete: {args.output_csv}")


if __name__ == "__main__":
    main()
