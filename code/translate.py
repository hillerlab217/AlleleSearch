# takes a fasta file and return two fasta files
# 3 parameters needed
# 1. input file path
# 2. output file path for longest orf
# 3. output file path for other orfs

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def translate_three_frames(dna_record):
    translations = {}
    for frame in range(3):
        subsequence = dna_record.seq[frame:]
        cleaned_sequence = subsequence.replace("-", "").replace("N", "")
        trimmed_sequence = cleaned_sequence[:len(cleaned_sequence) - (len(cleaned_sequence) % 3)]
        protein_sequence = trimmed_sequence.translate()
        translations[f"{dna_record.id}_frame{frame+1}"] = protein_sequence
    return translations

def extract_longest_orf(protein_seq):
    """Extracts the longest ORF starting with 'M' and ending with '*'"""
    orfs = []
    seq_str = str(protein_seq)
    start_index = -1

    for i, aa in enumerate(seq_str):
        if aa == "M" and start_index == -1:
            start_index = i
        if aa == "*" and start_index != -1:
            orfs.append(seq_str[start_index:i])
            start_index = -1

    if not orfs and start_index != -1:
        orfs.append(seq_str[start_index:])

    return max(orfs, key=len, default="")

def longest_reading_frame(translates):
    orfs = {frame: extract_longest_orf(seq) for frame, seq in translates.items()}
    sorted_orfs = sorted(orfs.items(), key=lambda x: len(x[1]), reverse=True)
    return {sorted_orfs[0][0]: translates[sorted_orfs[0][0]]}, {k: v for k, v in translates.items() if k != sorted_orfs[0][0]}

def main():
    parser = argparse.ArgumentParser(description="Translate DNA to protein and extract ORFs")
    parser.add_argument('input_file', help="Input DNA FASTA file")
    parser.add_argument('output_longest_orf', help="Output for longest ORF")
    parser.add_argument('output_others', help="Output for other ORFs")
    args = parser.parse_args()

    longest_readings = {}
    other_readings = {}

    for record in SeqIO.parse(args.input_file, "fasta"):
        translates = translate_three_frames(record)
        longest, others = longest_reading_frame(translates)
        # Check whether the best frame actually contains a valid ORF
        frame_id, frame_seq = next(iter(longest.items()))
        best_orf = extract_longest_orf(frame_seq)
        if not best_orf:
            print(
                f"WARNING: Skipping '{record.id}' — no valid ORF (Met→stop) "
                f"found in any reading frame. Sequence removed from output.",
                flush=True,
            )
            continue
        longest_readings.update(longest)
        other_readings.update(others)

    # Write outputs
    SeqIO.write((SeqRecord(Seq(v), id=k, description="") for k, v in longest_readings.items()), 
               args.output_longest_orf, "fasta")
    SeqIO.write((SeqRecord(Seq(v), id=k, description="") for k, v in other_readings.items()), 
               args.output_others, "fasta")

    print(f"Longest ORFs written to: {args.output_longest_orf}")
    print(f"Other ORFs written to: {args.output_others}")

if __name__ == "__main__":
    main()