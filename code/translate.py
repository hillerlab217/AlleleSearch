# takes a fasta file and return two fasta files
# 3 parameters needed
# 1. input file path
# 2. output file path for longest orf
# 3. output file path for other orfs

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def main():
    parser = argparse.ArgumentParser(description="Translate DNA to protein and extract ORFs")
    parser.add_argument('input_file', help="Input DNA FASTA file")
    parser.add_argument('output_longest_orf', help="Output for longest ORF")
    parser.add_argument('output_others', help="Output for other ORFs")
    args = parser.parse_args()

    longest_readings = {}
    other_readings = {}

    for record in SeqIO.parse(args.input_file, "fasta"):
        seq_str = str(record.seq).upper().replace('-', '').replace('N', '')
        seq_str = seq_str[:len(seq_str) - (len(seq_str) % 3)]
        if len(seq_str) < 3:
            print(
                f"WARNING: Skipping '{record.id}' — sequence too short to translate.",
                flush=True,
            )
            continue

        # codon_check upstream already validated and trimmed the DNA to begin at the
        # correct start codon, so Frame 0 is always the correct reading frame.
        # Replace the start codon with ATG so the protein always begins with M,
        # which is biologically correct (fMet initiation applies to all start codons).
        frame0_seq = 'ATG' + seq_str[3:]
        protein = str(Seq(frame0_seq).translate(to_stop=True))
        if not protein:
            print(
                f"WARNING: Skipping '{record.id}' — no protein sequence could be "
                f"extracted from Frame 0. Sequence removed from output.",
                flush=True,
            )
            continue
        longest_readings[record.id + '_frame1'] = protein

        # Frames 1 and 2 written to others for debugging; not used downstream.
        for frame in [1, 2]:
            s = seq_str[frame:]
            s = s[:len(s) - (len(s) % 3)]
            other_protein = str(Seq(s).translate(to_stop=True))
            if other_protein:
                other_readings[record.id + f'_frame{frame + 1}'] = other_protein

    # Write outputs
    SeqIO.write((SeqRecord(Seq(v), id=k, description="") for k, v in longest_readings.items()),
               args.output_longest_orf, "fasta")
    SeqIO.write((SeqRecord(Seq(v), id=k, description="") for k, v in other_readings.items()),
               args.output_others, "fasta")

    print(f"Longest ORFs written to: {args.output_longest_orf}")
    print(f"Other ORFs written to: {args.output_others}")

if __name__ == "__main__":
    main()