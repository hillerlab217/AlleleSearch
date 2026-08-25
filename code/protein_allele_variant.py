"""
Input: Takes two FASTA files (reference and strains), and a optional output directory path default as output/allele_info
Analysis:
    Identify allels, strains of each allels
    Identifies all variant positions
    Calculates number and percentage of variants
    Records amino acid changes
Output: CSV
Containing:
    Strain ID, number of strain, strain percentage
    Number of variants
    Percentage of variant positions
    Variant positions, reference, variant amino acids
"""
import sys
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner
from collections import defaultdict, Counter
import argparse
import csv
import re
from pathlib import Path
from header_utils import parse_header, species_from_id
from compute_consensus import compute_gapped_consensus_from_records

# allel_info list, for protein sequence
# allel_info = [str, int, float, int, float, defaultdict(list), list]
# [   0(str): allele name,
#     1(int): number of species have this allel,
#     2(float): percentage of species have this allel,
#     3(int): number of variant positions,
#     4(float): percentage of variant positions,
#     5(defaultdict): {variant position: [reference aa, variant aa]} position is 1-based,
#     6(list): [strain names that have this allele]
# ]


def allele_types(discovered_seqs, allele_suffix):
    """
    Analyze discovered sequences and create allele dictionary.
    Args:
        discovered_seqs: Dictionary of {sequence_id: protein_sequence}
    Returns:
        Dictionary {allele_sequence: allel_info_list}
    """
    # Group sequences by their protein sequence
    seq_to_strains = defaultdict(list)
    for strain_id, seq in discovered_seqs.items():
        seq_to_strains[seq].append(strain_id)
    
    total_sequences = len(discovered_seqs)
    
    # Sort alleles by count descending, then by sequence
    sorted_alleles = sorted(
        seq_to_strains.items(),
        key=lambda x: (-len(x[1]), x[0])
    )
    
    # Create allele dictionary
    allele_dict = {}
    for idx, (seq, strain_ids) in enumerate(sorted_alleles, 1):
        # Extract species names from the sequence ID part of each header
        species_names = [
            species_from_id(parse_header(sid)[0])
            for sid in strain_ids
        ]
        
        # Count unique species
        species_count = len(set(species_names))
        species_percentage = (species_count / total_sequences) * 100
        
        # Create allele info list
        allele_dict[seq] = [
            f"{allele_suffix}_{idx}",
            species_count,
            species_percentage,
            0,                             # placeholder
            0.0,                           # placeholder
            defaultdict(list),             # placeholder
            species_names
        ]
    
    return allele_dict


# Write to CSV
def write_csv(allele_dict, output_loc, with_strain_name=True, gapped_seqs=None, consensus_row=None):

    with open(output_loc, "w", newline="") as csvfile:
        fieldnames = [
            "allele_name", 
            "species_count", 
            "species_percentage", 
            "variant_position_count", 
            "variant_percentage",
            "variant_positions(position:reference→changed)", 
            "allele_sequence"
            ]
        
        if with_strain_name:
            fieldnames.append("strain_names")
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Prepend consensus row as first data row (no_strain CSV only, MAFFT mode)
        if consensus_row:
            label, gapped_seq = consensus_row
            writer.writerow({
                "allele_name": label,
                "species_count": "",
                "species_percentage": "",
                "variant_position_count": 0,
                "variant_percentage": "0.00",
                "variant_positions(position:reference→changed)": "",
                "allele_sequence": gapped_seq,
            })

        for allele_seq, allele_info in allele_dict.items():
            # Format variant positions dictionary as string
            variant_pos_str = "; ".join(
                f"{pos}:{vals[0]}→{vals[1]}"
                for pos, vals in sorted(allele_info[5].items())
            )

            # Use gapped sequence for the figure when available
            seq_for_csv = (gapped_seqs or {}).get(allele_seq, allele_seq)

            row_data = {
                "allele_name": allele_info[0],
                "species_count": allele_info[1],
                "species_percentage": f"{allele_info[2]:.2f}",
                "variant_position_count": allele_info[3],
                "variant_percentage": f"{allele_info[4]:.2f}",
                "variant_positions(position:reference→changed)":variant_pos_str,
                "allele_sequence": seq_for_csv
            }

            if with_strain_name:
                row_data["strain_names"] = "; ".join(sorted(allele_info[6]))

            writer.writerow(row_data)

def is_dna(sequence):
    aa_percentage = sum(1 for aa in sequence.upper() if aa in "ATCGatcg")/len(sequence)
    return aa_percentage > 0.8


# Find variation position
def variat_position(ref_seq, allele_dict):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    
    
    for seq, info in allele_dict.items():
        variants = defaultdict(list)
        ref_pos = 0
        num_diff = 0
        
        try:
            alignment = next(aligner.align(ref_seq, seq))
        except StopIteration:
            print(f"WARNING: Could not align allele '{info[0]}' to reference — skipping variant calling for this allele.", flush=True)
            info[3] = 0
            info[4] = 0.0
            info[5] = variants
            continue

        ref_aligned = alignment[0]
        seq_aligned = alignment[1]

        for ref_aa, seq_aa in zip(ref_aligned, seq_aligned):
            if ref_aa != '-':
                ref_pos += 1

            if ref_aa != seq_aa and seq_aa != '-':
                num_diff += 1
                variants[ref_pos] = [ref_aa, seq_aa]

        # Update allele info
        info[3] = num_diff
        info[4] = (num_diff / len(ref_seq)) * 100 if ref_seq else 0
        info[5] = variants
    
    return allele_dict


def variat_position_from_mafft(mafft_fasta, allele_dict):
    """Compute variant positions using an existing MAFFT multiple alignment.

    Positions are 1-based MAFFT alignment column numbers.
    Variant notation:
      substitution  →  pos : ref_aa  → allele_aa
      deletion      →  pos : ref_aa  → -          (allele has gap where consensus doesn't)
      insertion     →  pos : -       → allele_aa   (allele has AA where consensus has gap)

    The gapped allele sequence (with '-') is returned in a separate dict so
    the CSV can store it for the alignment figure.
    Returns (allele_dict, gapped_seqs, gapped_consensus_str).
    """
    try:
        records = list(AlignIO.read(mafft_fasta, "fasta"))
    except Exception as e:
        print(f"WARNING: Could not read MAFFT alignment '{mafft_fasta}': {e} — falling back to pairwise.", flush=True)
        return allele_dict, {}, ""

    if not records:
        return allele_dict, {}, ""

    alignment_len = len(records[0].seq)

    # Gapped consensus via shared helper (same algorithm as compute_consensus.py)
    gapped_consensus_str = compute_gapped_consensus_from_records(records)
    gapped_consensus = list(gapped_consensus_str)

    # Consensus length (non-gap positions) — used as denominator for variant %
    ref_len = sum(1 for aa in gapped_consensus if aa != '-')

    # Map ungapped protein sequence → gapped aligned sequence
    seq_to_gapped = {}
    for rec in records:
        gapped = str(rec.seq)
        ungapped = gapped.replace('-', '')
        seq_to_gapped[ungapped] = gapped

    gapped_seqs = {}   # ungapped_seq → gapped_seq, returned to caller

    for seq, info in allele_dict.items():
        variants = defaultdict(list)
        num_diff = 0

        gapped_seq = seq_to_gapped.get(seq)
        if gapped_seq is None:
            # Try stripping trailing stop codons
            stripped = seq.rstrip('*')
            for ug, g in seq_to_gapped.items():
                if ug.rstrip('*') == stripped:
                    gapped_seq = g
                    break

        if gapped_seq is None:
            print(f"WARNING: Allele '{info[0]}' not found in MAFFT alignment — skipping variant calling.", flush=True)
            info[3] = 0
            info[4] = 0.0
            info[5] = variants
            continue

        gapped_seqs[seq] = gapped_seq

        for col_idx in range(alignment_len):
            pos = col_idx + 1          # 1-based alignment column
            cons_aa = gapped_consensus[col_idx]
            allele_aa = gapped_seq[col_idx]

            if cons_aa == allele_aa:
                continue               # identical (including both '-')

            if cons_aa == '-' and allele_aa != '-':
                # Insertion in allele
                num_diff += 1
                variants[pos] = ['-', allele_aa]
            elif cons_aa != '-' and allele_aa == '-':
                # Deletion in allele
                num_diff += 1
                variants[pos] = [cons_aa, '-']
            else:
                # Substitution
                num_diff += 1
                variants[pos] = [cons_aa, allele_aa]

        info[3] = num_diff
        info[4] = (num_diff / ref_len) * 100 if ref_len else 0.0
        info[5] = variants

    return allele_dict, gapped_seqs, gapped_consensus_str


def most_frequent_allele(allele_dict):
    if not allele_dict:
        return None

    return max(allele_dict.items(), key=lambda x: x[1][1])[0]

# Write fasta file for allele

def write_allele_sequences_fasta(allele_dict, output_dir):
    output_path = Path(output_dir) / "allele_sequence_output.fasta"
    
    def extract_allele_number(info):
        match = re.search(r"_(\d+)$", info[0])
        return int(match.group(1)) if match else float('inf')

    # Sort allele_dict by allele number extracted from allele name
    sorted_alleles = sorted(
        allele_dict.items(),
        key=lambda x: extract_allele_number(x[1])
    )

    with open(output_path, "w") as fasta_out:
        for seq, info in sorted_alleles:
            header = f"{info[0]}|species_count={info[1]}|species_percentage={info[2]:.2f}"
            fasta_out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fasta_out.write(seq[i:i+60] + "\n")
    
    print(f"Allele sequences written to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Protein sequence variant analysis')
    parser.add_argument('-r', '--reference', default=None, help='Reference query FASTA file location, protein or dna sequence, '
                                            'if not given or None, will take the allele owned by most species as reference')
    parser.add_argument('discovered', help='Discovered protein FASTA file location')
    parser.add_argument('-o', '--output', default='output/allele_info', help='Output directory')
    parser.add_argument('-s', '--allele_suffix', default='allele', help='Allele name sufix')
    parser.add_argument('--mafft_alignment', default=None,
                        help='MAFFT-aligned FASTA of alleles; when provided, uses alignment '
                             'columns as positions and records insertions/deletions explicitly')
    args = parser.parse_args()

    from translate import translate_three_frames, longest_reading_frame

    # Create output directory if it doesn't exist
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    #get allele count
    try:
        discovered_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.discovered, 'fasta')}
        # Filter out zero-length sequences
        empty_ids = [sid for sid, seq in discovered_seqs.items() if len(seq) == 0]
        if empty_ids:
            print(f"WARNING: Skipping {len(empty_ids)} zero-length sequence(s): {empty_ids}", flush=True)
            discovered_seqs = {sid: seq for sid, seq in discovered_seqs.items() if len(seq) > 0}
        if not discovered_seqs:
            print("WARNING: No valid sequences to analyse — writing empty output files.", flush=True)
            output_dir = Path(args.output)
            output_dir.mkdir(parents=True, exist_ok=True)
            for fname in ("allele_info_with_strain_name.csv", "allele_info_no_strain_name.csv"):
                with open(output_dir / fname, "w", newline="") as f:
                    f.write("allele_name,species_count,species_percentage,variant_position_count,variant_percentage,variant_positions(position:reference→changed),allele_sequence\n")
            return
        allele_dict = allele_types(discovered_seqs, args.allele_suffix)
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)

    # Read sequences
    effective_mafft = args.mafft_alignment  # may be set to None in empty-FASTA fallback
    if args.reference is not None:
        recs = list(SeqIO.parse(args.reference, 'fasta'))
        if recs:
            ref_rec = recs[0]
            ref_seq = str(ref_rec.seq)
            if is_dna(ref_seq):
                translates = translate_three_frames(ref_rec)
                ref_seq, _ = longest_reading_frame(translates)
                ref_seq = str(list(ref_seq.values())[0])
        else:
            print(
                f"WARNING: Consensus FASTA '{args.reference}' is empty. "
                f"Falling back to self-reference (most frequent allele). "
                f"Re-run the full pipeline to regenerate the consensus.",
                file=sys.stderr,
            )
            ref_seq = most_frequent_allele(allele_dict)
            effective_mafft = None  # force pairwise mode
    else:
        ref_seq = most_frequent_allele(allele_dict)

    try:
        if effective_mafft:
            allele_dict, gapped_seqs, gapped_consensus_str = variat_position_from_mafft(effective_mafft, allele_dict)
        else:
            allele_dict = variat_position(ref_seq, allele_dict)
            gapped_seqs = {}
            gapped_consensus_str = None

        # Write the fasta file for alleles (always ungapped sequences)
        write_allele_sequences_fasta(allele_dict, output_dir)

        # Generate output file paths
        output_with_strain = output_dir / "allele_info_with_strain_name.csv"
        output_no_strain = output_dir / "allele_info_no_strain_name.csv"

        # Embed gapped consensus as first row in no-strain CSV (consensus_reference mode only)
        consensus_label = f"{args.allele_suffix}_consensus" if gapped_consensus_str else None
        write_csv(allele_dict, output_with_strain, with_strain_name=True, gapped_seqs=gapped_seqs)
        write_csv(allele_dict, output_no_strain, with_strain_name=False, gapped_seqs=gapped_seqs,
                  consensus_row=(consensus_label, gapped_consensus_str) if consensus_label else None)

        print(f"Analysis complete. Results saved to: {output_dir}")

    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()