"""
extend_to_stop.py

For each sequence in a no_stop FASTA (has start codon, missing stop codon),
scans the genome downstream in-frame until a stop codon is found.
Outputs a FASTA with the same header format as standard_dna_seq:

  >ID|contig=...|E-value=...|AlignLen=...|SeqIdentity=...|SubjectStart=...|SubjectEnd=...|StartCodonOffset=...|StopCodonOffset=...

Usage:
    python extend_to_stop.py no_stop.fasta /path/to/genome_dir output.fasta
    python extend_to_stop.py no_stop.fasta /path/to/genome_dir output.fasta --max_codons 200
"""

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from header_utils import parse_header, build_header, species_from_id

STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def get_downstream_extension(genome_filepath, contig_name, subject_start, subject_end, max_extend_bp):
    """
    Fetch the downstream genome sequence beyond SubjectEnd, in gene orientation (5'→3').

    Forward strand (subject_start < subject_end): extension is bases after subject_end.
    Reverse strand (subject_start > subject_end): extension is bases before subject_end,
    reverse-complemented.

    Returns (extension_seq, is_reverse) or (None, None) on failure.
    """
    if not os.path.exists(genome_filepath):
        print(f"Warning: genome file not found: {genome_filepath}", file=sys.stderr)
        return None, None

    with open(genome_filepath) as f:
        for record in SeqIO.parse(f, "fasta"):
            if contig_name != record.id and contig_name not in record.id:
                continue

            genome = str(record.seq).upper()
            genome_len = len(genome)
            is_reverse = subject_start > subject_end

            if is_reverse:
                # Downstream in gene direction = lower genome coordinates.
                # subject_end is 1-based; the extension occupies coords below it.
                ext_end_0   = subject_end - 1
                ext_start_0 = max(0, ext_end_0 - max_extend_bp)
                fragment    = genome[ext_start_0:ext_end_0]
                extension   = str(Seq(fragment).reverse_complement())
            else:
                # Downstream = higher genome coordinates beyond subject_end (1-based inclusive).
                ext_start_0 = subject_end
                ext_end_0   = min(genome_len, subject_end + max_extend_bp)
                extension   = genome[ext_start_0:ext_end_0]

            return extension, is_reverse

    print(f"Warning: contig '{contig_name}' not found in {genome_filepath}", file=sys.stderr)
    return None, None


def scan_for_stop(seq, extension):
    """
    Concatenate seq + extension and scan in-frame (step 3 from position 0) for the
    first stop codon.  seq is guaranteed to have no stop codons (no_stop classification).

    Returns (trimmed_seq, bases_added) or (None, None) if no stop found.
    """
    full = seq + extension
    for i in range(0, len(full) - 2, 3):
        if full[i:i+3] in STOP_CODONS:
            trimmed = full[:i + 3]
            return trimmed, len(trimmed) - len(seq)
    return None, None


def main():
    parser = argparse.ArgumentParser(
        description="Extend no_stop sequences downstream until a stop codon is found."
    )
    parser.add_argument("input_fasta",  help="Input FASTA file (no_stop sequences from codon check)")
    parser.add_argument("genome_dir",   help="Directory containing genome FASTA files")
    parser.add_argument("output_fasta", help="Output FASTA file")
    parser.add_argument("--max_codons", type=int, default=500,
                        help="Maximum codons to scan downstream (default: 500)")
    args = parser.parse_args()

    max_extend_bp = args.max_codons * 3
    written = skipped = 0

    with open(args.output_fasta, "w") as out_f:
        with open(args.input_fasta) as in_f:
            for record in SeqIO.parse(in_f, "fasta"):
                header = ">" + record.description
                seq    = str(record.seq).upper()
                seq_id, metadata = parse_header(header)

                missing = [f for f in ('contig', 'SubjectStart', 'SubjectEnd') if f not in metadata]
                if missing:
                    print(f"Skipping {seq_id}: missing fields {missing}", file=sys.stderr)
                    skipped += 1
                    continue

                contig        = metadata['contig']
                subject_start = int(metadata['SubjectStart'])
                subject_end   = int(metadata['SubjectEnd'])
                species       = species_from_id(seq_id)
                # Genome files may have version/assembly suffixes (e.g. GCA_001234.2_ASM_genomic.fasta).
                # Find any file whose name starts with the accession prefix.
                matches = [f for f in os.listdir(args.genome_dir)
                           if f.startswith(species) and f.endswith(".fasta")]
                if not matches:
                    print(f"Skipping {seq_id}: no genome file starting with '{species}' in {args.genome_dir}",
                          file=sys.stderr)
                    skipped += 1
                    continue
                genome_file = os.path.join(args.genome_dir, matches[0])

                extension, is_reverse = get_downstream_extension(
                    genome_file, contig, subject_start, subject_end, max_extend_bp
                )
                if extension is None:
                    skipped += 1
                    continue

                trimmed, bases_added = scan_for_stop(seq, extension)
                if trimmed is None:
                    print(f"Skipping {seq_id}: no stop codon within {args.max_codons} codons downstream",
                          file=sys.stderr)
                    skipped += 1
                    continue

                # Update header to match standard_dna_seq format
                metadata['AlignLen'] = str(len(trimmed))
                if is_reverse:
                    metadata['SubjectEnd'] = str(subject_end - bases_added)
                else:
                    metadata['SubjectEnd'] = str(subject_end + bases_added)
                # StopCodonOffset: codons extended beyond original BLAST end
                metadata['StopCodonOffset'] = str(bases_added // 3)
                new_header = build_header(seq_id, metadata)

                out_f.write(f"{new_header}\n")
                for i in range(0, len(trimmed), 60):
                    out_f.write(f"{trimmed[i:i+60]}\n")
                written += 1

    print(f"Done. Written: {written}, skipped: {skipped}")


if __name__ == "__main__":
    main()
