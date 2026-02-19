"""
look for start codon and stop codon infame, 
input a dna sequence fasta file, a output directory name 
return dna fasta files according to the type of dna sequence

change the start and stop site according to the trimm

    "standard" have start and stop codon in frame
    "no_stop_no_start"
    "no_start"
    "no_stop"
    "no_start_stop_not_at_end"
    "stop_not_at_end"
    "more_than_one_stop"

 """

import os
import argparse
from collections import OrderedDict
from Bio.Seq import Seq
from header_utils import parse_header, build_header

# Define codon sets
start_codons = {'ATG', 'GTG', 'TTG'}
stop_codons = {'TAA', 'TAG', 'TGA'}

def have_start_codon(dna_seq, codon_range=10, direction="down"):
    """Check for start codon in the first `codon_range` codons."""
    if direction=="down":
        for i in range(min(codon_range, len(dna_seq) // 3)):
            codon = dna_seq[i*3 : (i+1)*3]
            if codon in start_codons:
                return True, i*3  # Return position of the start codon
        return False, 0
    elif direction=="up":
        for i in reversed(range(min(codon_range, len(dna_seq) // 3))):
            codon = dna_seq[i*3 : (i+1)*3]
            if codon in start_codons:
                return True, i*3  # Return position of the start codon
        return False, 0
    else:
        raise ValueError("direction must be either 'down' or 'up'")


def end_with_stop_codon(dna_seq, codon_range=10):
    """Check for stop codon in the last `codon_range` codons."""
    seq_len = len(dna_seq)
    left_over = seq_len%3
    for i in reversed(range(1, min(codon_range, seq_len // 3) + 1)):
        codon = dna_seq[seq_len - left_over - i*3 : seq_len -left_over- (i-1)*3]
        if codon in stop_codons:
            return True, seq_len - (i-1)*3  # Return end position of the stop codon
    return False, seq_len

def count_stop_codons(dna_seq, codon_range=10):
    """Count in-frame stop codons (excluding trailing partial codons)."""
    seq_len = len(dna_seq)
    count = 0
    # Search the entire sequence except the last codon_range
    search_end = seq_len - (seq_len % 3) - (codon_range*3 if codon_range*3 < seq_len else 0)
    for i in range(0, max(0, search_end), 3):
        if dna_seq[i:i+3] in stop_codons:
            count += 1
    return count

def classify_sequence(dna_seq, codon_range=10):
    """Classify sequence based on start/stop codon presence."""
    has_start, start_pos = have_start_codon(dna_seq, codon_range)
    has_stop, stop_pos = end_with_stop_codon(dna_seq, codon_range)
    stop_count = count_stop_codons(dna_seq, codon_range)

    if stop_count > 1:
        return "more_than_one_stop", start_pos, stop_pos
    
    if has_start and has_stop and stop_count == 0:
        return "standard", start_pos, stop_pos
    
    if stop_count == 1 and not has_stop:
        if has_start:
            return "stop_not_at_end", start_pos, len(dna_seq)
        else:
            return "no_start_stop_not_at_end", 0, len(dna_seq)
    
    if not has_stop and stop_count == 0:
        if has_start:
            return "no_stop", start_pos, len(dna_seq)
        else:
            return "no_stop_no_start", 0, len(dna_seq)
    
    if has_stop and not has_start:
        return "no_start", 0, stop_pos
    
    # Default case (shouldn't normally reach here)
    return "no_stop_no_start", 0, len(dna_seq)

def check_stop_start_codon_dna(dna_string, codon_range=10):
    dna_string = dna_string.replace('-', '')

    status_que = {
        "standard": 1,
        "no_stop_no_start": 2,
        "no_stop": 3,
        "no_start": 4,
        "no_start_stop_not_at_end": 5,
        "stop_not_at_end": 6,
        "more_than_one_stop": 7
    }

    status = "more_than_one_stop"
    final_start = 0
    final_stop = len(dna_string)

    for i in range(3):
        dna_seq = dna_string[i:]
        status_temp, start_position, stop_position = classify_sequence(dna_seq, codon_range)
        if status_que[status] > status_que[status_temp]:
            status = status_temp
            final_start = start_position + i
            final_stop = stop_position + i

    return status, final_start, final_stop

def update_header(header, new_start, new_len, is_reverse=False):
    """Update subject coordinates and alignment length after sequence trimming.

    new_start : offset (in bases) from the original SubjectStart where the
                trimmed sequence begins.
    new_len   : length of the trimmed sequence in bases.
    is_reverse: True when the hit is on the reverse strand (SubjectStart >
                SubjectEnd in the original BLAST result).
    """
    seq_id, metadata = parse_header(header)

    if 'SubjectStart' in metadata and 'SubjectEnd' in metadata:
        try:
            orig_start = int(metadata['SubjectStart'])
            if is_reverse:
                new_subj_start = orig_start - new_start
                new_subj_end = new_subj_start - new_len + 1
            else:
                new_subj_start = orig_start + new_start
                new_subj_end = new_subj_start + new_len - 1
            metadata['SubjectStart'] = str(new_subj_start)
            metadata['SubjectEnd'] = str(new_subj_end)
        except ValueError:
            pass

    metadata['AlignLen'] = str(new_len)
    return build_header(seq_id, metadata)

def process_fasta(input_file, output_dir, codon_range=10):
    """Classify and trim sequences from input FASTA."""
    os.makedirs(output_dir, exist_ok=True)
    sequences = {status: [] for status in [
        "standard", "no_stop_no_start", "no_start", "no_stop",
        "no_start_stop_not_at_end", "stop_not_at_end", "more_than_one_stop"
    ]}

    with open(input_file) as f:
        current_id, current_seq = None, []
        for line in f:
            if line.startswith('>'):
                if current_id:
                    process_sequence(current_id, ''.join(current_seq), sequences, codon_range)
                current_id = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            process_sequence(current_id, ''.join(current_seq), sequences, codon_range)

    # Write output files
    for status, seqs in sequences.items():
        if seqs:
            with open(f"{output_dir}/{status}.fasta", 'w') as f:
                for header, seq in seqs:
                    f.write(f"{header}\n")
                    for i in range(0, len(seq), 60):
                        f.write(f"{seq[i:i+60]}\n")

def process_sequence(header, dna_seq, sequences, codon_range):
    """Classify and trim a single sequence."""
    dna_seq = dna_seq.upper().replace('-', '')
    status, start_pos, end_pos = check_stop_start_codon_dna(dna_seq, codon_range)
    
    # Determine strand direction from header
    is_reverse = False
    _, metadata = parse_header(header)
    if 'SubjectStart' in metadata and 'SubjectEnd' in metadata:
        try:
            subj_start = int(metadata['SubjectStart'])
            subj_end = int(metadata['SubjectEnd'])
            is_reverse = subj_start > subj_end
        except ValueError:
            pass
    
    # Trim sequence
    trimmed_seq = dna_seq[start_pos:end_pos]

    
    # Update header
    new_header = update_header(
        header, start_pos, len(trimmed_seq), is_reverse
    )
    sequences[status].append((new_header, trimmed_seq))

def main():
    parser = argparse.ArgumentParser(
        description="Classify DNA sequences by start/stop codon presence."
    )
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--codon_range", type=int, default=10,
                       help="Number of codons to check for start/stop (default: 10)")
    parser.add_argument("--start_search_direction", type=str, default="down",
                       help="direction must be either 'down' or 'up'(default: down)")
    args = parser.parse_args()
    process_fasta(args.input_fasta, args.output_dir, args.codon_range)

if __name__ == "__main__":
    main()