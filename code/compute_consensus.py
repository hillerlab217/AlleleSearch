"""
Compute a consensus sequence from a MAFFT-aligned FASTA file.
Uses BioPython's SummaryInfo.gap_consensus: at each alignment column picks
the most common character (including gaps) if it appears in >=50% of all
sequences; gaps count in both the vote and the denominator, so a position
where the majority of sequences have '-' yields '-' in the consensus.
Ambiguous positions (no character reaches 50%) produce 'X'.
All-gap columns produce '-'.
"""
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
import argparse


def compute_gapped_consensus_from_records(records):
    """Return a gapped consensus string of length == alignment_len using
    BioPython's gap_consensus (threshold=0.5, gaps included in vote).
    Positions where no character reaches 50% produce 'X'.
    All-gap columns produce '-'. Importable by other scripts."""
    alignment = MultipleSeqAlignment(records)
    summary = AlignInfo.SummaryInfo(alignment)
    consensus = str(summary.gap_consensus(threshold=0.5, ambiguous='X'))
    # Replace 'X' with '-' for purely all-gap columns (already '-' via
    # gap_consensus, but kept as a safety net)
    result = list(consensus)
    for i in range(alignment.get_alignment_length()):
        if all(str(rec.seq[i]) == '-' for rec in alignment):
            result[i] = '-'
    return ''.join(result)


def compute_consensus(aligned_fasta, name, output_fasta):
    try:
        records = list(AlignIO.read(aligned_fasta, "fasta"))
    except (ValueError, FileNotFoundError):
        records = []

    if len(records) == 0:
        open(output_fasta, 'w').close()
        print(f"No sequences in alignment; writing empty consensus file: {output_fasta}", flush=True)
        return

    if len(records) == 1:
        # Only one allele — it is its own consensus
        seq = str(records[0].seq).replace('-', '')
    else:
        gapped = compute_gapped_consensus_from_records(records)
        seq = gapped.replace('-', '')

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
