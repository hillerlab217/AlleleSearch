"""
Compute a consensus sequence from a MAFFT-aligned FASTA file.
At each alignment column, picks the most common residue (ignoring gaps).
"""
from Bio import AlignIO
from collections import Counter
import argparse


def compute_consensus(aligned_fasta, name, output_fasta):
    records = list(AlignIO.read(aligned_fasta, "fasta"))

    if len(records) == 0:
        raise ValueError("No sequences found in alignment file")

    if len(records) == 1:
        # Only one allele — it is its own consensus
        seq = str(records[0].seq).replace('-', '')
    else:
        alignment_len = len(records[0].seq)
        consensus = []
        for i in range(alignment_len):
            col = [str(rec.seq[i]) for rec in records]
            non_gap = [aa for aa in col if aa != '-']
            if non_gap:
                consensus.append(Counter(non_gap).most_common(1)[0][0])
        seq = ''.join(consensus)

    with open(output_fasta, 'w') as f:
        f.write(f">{name}_consensus\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + '\n')

    print(f"Consensus sequence ({len(seq)} aa) written to: {output_fasta}", flush=True)


def main():
    parser = argparse.ArgumentParser(description='Compute consensus from MAFFT alignment')
    parser.add_argument('aligned_fasta', help='MAFFT-aligned FASTA file')
    parser.add_argument('-n', '--name', default='consensus', help='Prefix for the consensus FASTA header')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file path')
    args = parser.parse_args()

    compute_consensus(args.aligned_fasta, args.name, args.output)


if __name__ == '__main__':
    main()
