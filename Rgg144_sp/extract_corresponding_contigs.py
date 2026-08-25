"""
Extract the contig containing Rgg144 from each species genome in GoldenSet.

Inputs:
  - standard_dna_seq/Rgg144.fasta       : FASTA with headers encoding strain → contig name
  - allele_info/.../allele_info_with_strain_name.csv : strain names per allele
  - /home/xinruoz/Document/Hiller/data/GoldenSet : genome FASTA files

Output:
  - corresponding_contigs/{strain}.fasta : one file per species, containing the
    contig that holds the Rgg144 gene
"""

import argparse
import csv
import os
import sys
from pathlib import Path

from Bio import SeqIO


def parse_strain_contig_map(standard_fasta):
    """Parse standard_dna_seq FASTA headers -> dict[strain] = contig_name."""
    strain_contig = {}
    with open(standard_fasta) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            header = line[1:].strip()
            # Format: Rgg144_vs_{strain}|contig={contig}|...
            parts = header.split('|')
            id_field = parts[0]  # e.g. Rgg144_vs_1235_3063_00
            # strain is everything after first 'vs_'
            vs_idx = id_field.find('vs_')
            if vs_idx == -1:
                continue
            strain = id_field[vs_idx + 3:]
            # contig field
            contig = None
            for part in parts[1:]:
                if part.startswith('contig='):
                    contig = part[len('contig='):]
                    break
            if contig:
                strain_contig[strain] = contig
    return strain_contig


def parse_strains_from_allele_csv(csv_path):
    """Return sorted list of unique strain names from allele_info CSV."""
    strains = set()
    with open(csv_path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            names = row.get('strain_names', '')
            for name in names.split(';'):
                name = name.strip()
                if name:
                    strains.add(name)
    return sorted(strains)


def contig_matches(record_id, contig_name):
    """True if the genome contig record matches contig_name prefix-safely."""
    if record_id == contig_name:
        return True
    if record_id.startswith(contig_name + '|'):
        return True
    if record_id.startswith(contig_name + ' '):
        return True
    return False


def extract_contig(genome_fasta, contig_name):
    """Return list of SeqRecord whose ID matches contig_name."""
    matches = []
    for record in SeqIO.parse(genome_fasta, 'fasta'):
        if contig_matches(record.id, contig_name):
            matches.append(record)
    return matches


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--standard_fasta',
        default='standard_dna_seq/Rgg144.fasta',
        help='standard_dna_seq FASTA with strain/contig info in headers',
    )
    parser.add_argument(
        '--allele_csv',
        default='allele_info/Rgg144/self_reference/allele_info_with_strain_name.csv',
        help='allele_info CSV with strain_names column',
    )
    parser.add_argument(
        '--genome_dir',
        default='/home/xinruoz/Document/Hiller/data/GoldenSet',
        help='Directory containing genome FASTA files ({strain}.fasta)',
    )
    parser.add_argument(
        '--outdir',
        default='corresponding_contigs',
        help='Output directory (one FASTA per species)',
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f'Parsing contig info from: {args.standard_fasta}')
    strain_contig = parse_strain_contig_map(args.standard_fasta)
    print(f'  Found {len(strain_contig)} strains with contig info')

    print(f'Parsing strain names from: {args.allele_csv}')
    strains = parse_strains_from_allele_csv(args.allele_csv)
    print(f'  Found {len(strains)} unique strains in allele info')

    missing_contig = []
    missing_genome = []
    no_contig_found = []
    written = 0

    for strain in strains:
        # Get contig name
        if strain not in strain_contig:
            missing_contig.append(strain)
            continue

        contig_name = strain_contig[strain]

        # Find genome FASTA
        genome_fasta = Path(args.genome_dir) / f'{strain}.fasta'
        if not genome_fasta.exists():
            missing_genome.append(strain)
            continue

        # Extract contig
        records = extract_contig(str(genome_fasta), contig_name)
        if not records:
            no_contig_found.append((strain, contig_name))
            continue

        # Write output
        out_path = Path(args.outdir) / f'{strain}.fasta'
        SeqIO.write(records, str(out_path), 'fasta')
        written += 1

    print(f'\nDone. Written: {written} files -> {args.outdir}/')
    if missing_contig:
        print(f'  Skipped (no contig info in standard_fasta): {len(missing_contig)}')
        for s in missing_contig[:10]:
            print(f'    {s}')
        if len(missing_contig) > 10:
            print(f'    ... and {len(missing_contig) - 10} more')
    if missing_genome:
        print(f'  Skipped (genome FASTA not found in GoldenSet): {len(missing_genome)}')
        for s in missing_genome[:10]:
            print(f'    {s}')
        if len(missing_genome) > 10:
            print(f'    ... and {len(missing_genome) - 10} more')
    if no_contig_found:
        print(f'  Skipped (contig not found in genome file): {len(no_contig_found)}')
        for s, c in no_contig_found[:10]:
            print(f'    {s} (contig={c})')
        if len(no_contig_found) > 10:
            print(f'    ... and {len(no_contig_found) - 10} more')


if __name__ == '__main__':
    main()
