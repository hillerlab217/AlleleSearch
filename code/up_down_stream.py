import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from header_utils import parse_header, build_header, species_from_id

def find_extended_sequence(mode, genome_index, contig_name, start_position, end_position, downstream_length):
    """
    Return the extended sequence from the specified contig_name within the
    pre-loaded genome_index (dict[contig_id, SeqRecord]), from start_position
    to end_position plus the downstream_length.
    """

    genome_record = genome_index.get(contig_name)
    if genome_record is None:
        # Fallback: partial match (contig_name contained in record ID)
        matches = [r for rid, r in genome_index.items() if contig_name in rid]
        genome_record = matches[0] if matches else None
    if genome_record is None:
        print(f"Warning: Contig '{contig_name}' not found in genome index.")
        return None, None, None

    genome_sequence = str(genome_record.seq)
    genome_length = len(genome_sequence)

    if start_position < end_position:  # Forward strand
        if mode == 'down_stream':
            seq_start = start_position - 1  # 0 based
            seq_end = min(end_position + downstream_length, genome_length)
        elif mode == 'up_stream':
            seq_start = max(start_position - 1 - downstream_length, 0)
            seq_end = end_position
        elif mode == 'up_down_stream':
            seq_start = max(start_position - 1 - downstream_length, 0)
            seq_end = min(end_position + downstream_length, genome_length)
        else:
            raise ValueError("Mode not included, choose mode of program: down_stream, up_stream, up_down_stream")

        extended_seq = genome_sequence[seq_start:seq_end]
        seq_start = seq_start + 1

    else:  # Reverse strand
        if mode == 'down_stream':
            seq_start = start_position  # 0 based
            seq_end = max(0, end_position - 1 - downstream_length)
        elif mode == 'up_stream':
            seq_start = min(genome_length, start_position + downstream_length)
            seq_end = end_position - 1
        elif mode == 'up_down_stream':
            seq_start = min(genome_length, start_position + downstream_length)
            seq_end = max(0, end_position - 1 - downstream_length)
        else:
            raise ValueError("Mode not included, choose mode of program: down_stream, up_stream, up_down_stream")

        fragment = genome_sequence[seq_end:seq_start]
        extended_seq = str(Seq(fragment).reverse_complement())
        seq_end = seq_end + 1

    return extended_seq, seq_start, seq_end
def format_fasta_sequence(sequence, line_length=60):
    """Formats a sequence so that it wraps every `line_length` characters."""
    return "\n".join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))


def update_header(header, new_len, new_start, new_end, mode, extension_length):
    """Update a header with absolute subject coordinates and a new alignment length.

    Unlike the version in start_stop_codon_dna_seq.py, new_start and new_end
    are the final absolute genome positions (not offsets), computed by the
    caller after extending the sequence.

    Also records UpstreamExtension and DownstreamExtension (in base pairs) so
    that start_stop_codon_dna_seq.py can compute codon offsets relative to the
    original BLAST hit boundary.
    """
    seq_id, metadata = parse_header(header)
    if 'SubjectStart' in metadata and 'SubjectEnd' in metadata:
        metadata['SubjectStart'] = str(new_start)
        metadata['SubjectEnd'] = str(new_end)
    metadata['AlignLen'] = str(new_len)
    if mode == 'up_stream':
        metadata['UpstreamExtension'] = str(extension_length)
        metadata['DownstreamExtension'] = '0'
    elif mode == 'down_stream':
        metadata['UpstreamExtension'] = '0'
        metadata['DownstreamExtension'] = str(extension_length)
    elif mode == 'up_down_stream':
        metadata['UpstreamExtension'] = str(extension_length)
        metadata['DownstreamExtension'] = str(extension_length)
    return build_header(seq_id, metadata)

def main():
    parser = argparse.ArgumentParser(description="Find extended sequences for protein sequences.")
    parser.add_argument("--mode", required=True, help="Choose mode: down_stream, up_stream, up_down_stream")
    parser.add_argument("--seq_file", required=True, help="Path to the protein FASTA file.")
    parser.add_argument("--genome_dir", required=True, help="Directory containing genome FASTA files.")
    parser.add_argument("--length", type=int, required=True, help="Length of upstream/downstream sequence to extract.")
    parser.add_argument("--min_aa_length", type=int, default=5, help="Minimum amino acid length to keep a sequence (default: 5).")
    parser.add_argument("--output_fasta", type=str, required=True, help="Output FASTA file path.")
    args = parser.parse_args()

    mode = args.mode
    seq_filepath = args.seq_file
    genome_dir = args.genome_dir
    downstream_length = args.length
    min_aa_length = args.min_aa_length
    output_filepath = args.output_fasta

    strains = {}

    with open(seq_filepath, "r") as protein_handle:
        for protein_record in SeqIO.parse(protein_handle, "fasta"):
            header = protein_record.description.strip()
            strain_name, metadata = parse_header(header)

            required_fields = ['contig', 'E-value', 'AlignLen', 'SeqIdentity', 'SubjectStart', 'SubjectEnd']
            missing_fields = [f for f in required_fields if f not in metadata]

            if missing_fields:
                print(f"Skipping header (missing fields: {missing_fields}): {header}")
                continue

            if strain_name not in strains:
                strains[strain_name] = []

            strains[strain_name].append({
                "header": header,
                "contig": metadata['contig'],
                "e_value": metadata['E-value'],
                "align_len": metadata['AlignLen'],
                "seq_identity": metadata['SeqIdentity'],
                "subject_start": int(metadata['SubjectStart']),
                "subject_end": int(metadata['SubjectEnd'])
            })

    # Process each strain
    for strain_name, entries in strains.items():
        actual_strain_name = species_from_id(strain_name)
        # Genome files may have version/assembly suffixes (e.g. GCA_001234.2_ASM_genomic.fasta).
        # Find any file whose name starts with the accession prefix.
        matches = [f for f in os.listdir(genome_dir)
                   if f.startswith(actual_strain_name) and f.endswith(".fasta")]
        if not matches:
            print(f"Warning: No genome file starting with '{actual_strain_name}' in {genome_dir}. Skipping.")
            continue
        genome_filepath = os.path.join(genome_dir, matches[0])

        # Parse genome once per strain and index by contig ID
        genome_index = SeqIO.to_dict(SeqIO.parse(genome_filepath, "fasta"))

        for entry in entries:
            if entry["subject_start"] < 0 or entry["subject_end"] < 0:
                print(f"Warning: Negative positions for {strain_name}. Skipping.")
                continue

            result = find_extended_sequence(
                mode,
                genome_index,
                entry["contig"],
                entry["subject_start"],
                entry["subject_end"],
                downstream_length
            )

            if not result:
                print(f"Warning: Sequence not found for contig {entry['contig']} in {genome_filepath}")
                continue

            extended_seq, new_start, new_end = result

            if len(extended_seq) < min_aa_length * 3:
                print(f"Skipping sequence (length {len(extended_seq)} bp < {min_aa_length} aa minimum): {entry['header']}")
                continue

            entry["extended_seq"] = extended_seq

            updated_header = update_header(
                entry["header"],
                new_len=abs(new_end - new_start)+1,
                new_start=new_start,
                new_end=new_end,
                mode=mode,
                extension_length=downstream_length,
            )
            entry["updated_header"] = updated_header

    # Write output FASTA
    with open(output_filepath, "w") as output_handle:
        for strain_name, entries in strains.items():
            for entry in entries:
                if "extended_seq" in entry:
                    formatted_seq = format_fasta_sequence(entry['extended_seq'])
                    output_handle.write(f"{entry['updated_header']}\n{formatted_seq}\n")

    print(f"Extended sequences written to {output_filepath}")

if __name__ == "__main__":
    main()