import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from header_utils import parse_header, build_header, species_from_id


def extract_region(genome_filepath, contig, start, end, upstream_bp, downstream_bp):
    """
    Extract upstream and downstream flanking regions around a gene hit.

    start, end: 1-based BLAST coordinates. start > end signals reverse strand.
    upstream_bp / downstream_bp: bp to extract in the gene's 5' / 3' direction.

    Returns (upstream_seq, gene_seq, downstream_seq) all oriented 5'→3' relative
    to the gene (reverse-complemented when the hit is on the reverse strand).
    Any region that reaches a genome boundary is silently clamped.
    Returns (None, None, None) if the genome file or contig is not found.
    """
    if not os.path.exists(genome_filepath):
        print(f"Warning: Genome file {genome_filepath} not found. Skipping.")
        return None, None, None

    forward = start < end
    if forward:
        gene_start0 = start - 1  # convert to 0-based
        gene_end0 = end           # exclusive upper bound
    else:
        gene_start0 = end - 1    # smaller coordinate (0-based)
        gene_end0 = start         # exclusive upper bound

    with open(genome_filepath, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            if contig == record.id or contig in record.id:
                seq = str(record.seq)
                glen = len(seq)

                up_raw   = seq[max(0, gene_start0 - upstream_bp) : gene_start0]
                gene_raw = seq[gene_start0 : gene_end0]
                dn_raw   = seq[gene_end0   : min(glen, gene_end0 + downstream_bp)]

                if forward:
                    return up_raw, gene_raw, dn_raw
                else:
                    # Reverse strand: swap upstream/downstream and reverse-complement
                    return (
                        str(Seq(dn_raw).reverse_complement()),
                        str(Seq(gene_raw).reverse_complement()),
                        str(Seq(up_raw).reverse_complement()),
                    )

    print(f"Warning: Contig '{contig}' not found in {genome_filepath}.")
    return None, None, None


def format_fasta(sequence, line_length=60):
    return "\n".join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))


def write_record(fh, header, seq_id, metadata, region_label, sequence):
    metadata["Region"] = region_label
    new_header = build_header(seq_id, metadata)
    fh.write(f"{new_header}\n{format_fasta(sequence)}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract upstream/downstream flanking regions for gene hits."
    )
    parser.add_argument("--seq_file",    required=True,  help="Input FASTA (standard_dna_seq format).")
    parser.add_argument("--genome_dir",  required=True,  help="Directory of genome FASTA files.")
    parser.add_argument("--upstream",    type=int, default=0, help="bp upstream to extract (default: 0).")
    parser.add_argument("--downstream",  type=int, default=0, help="bp downstream to extract (default: 0).")
    parser.add_argument("--flank_only",  action="store_true",  help="Exclude the gene itself from output.")
    parser.add_argument("--merge",       action="store_true",  help="Merge upstream+gene+downstream into one file.")
    parser.add_argument("--output_dir",  default=".",    help="Directory to write output files (default: cwd).")
    args = parser.parse_args()

    upstream_bp   = args.upstream
    downstream_bp = args.downstream
    flank_only    = args.flank_only
    merge         = args.merge
    output_dir    = args.output_dir

    if upstream_bp == 0 and downstream_bp == 0:
        print("Warning: Both --upstream and --downstream are 0. Nothing to extract.")
        sys.exit(0)

    # Warn if merge requested but only one direction specified
    if merge and (upstream_bp == 0 or downstream_bp == 0):
        present = "upstream" if upstream_bp > 0 else "downstream"
        print(f"Warning: --merge requested but only --{present} is set. Falling back to split output.")
        merge = False

    os.makedirs(output_dir, exist_ok=True)

    # Open output file handles
    up_fh     = open(os.path.join(output_dir, "upstream.fasta"),   "w") if upstream_bp   > 0 and not merge else None
    dn_fh     = open(os.path.join(output_dir, "downstream.fasta"), "w") if downstream_bp > 0 and not merge else None
    merged_fh = open(os.path.join(output_dir, "merged.fasta"),     "w") if merge else None

    skipped = 0

    with open(args.seq_file, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            header = record.description.strip()
            seq_id, metadata = parse_header(header)

            required = ["contig", "SubjectStart", "SubjectEnd"]
            missing  = [f for f in required if f not in metadata]
            if missing:
                print(f"Skipping record (missing fields {missing}): {header}")
                skipped += 1
                continue

            contig = metadata["contig"]
            start  = int(metadata["SubjectStart"])
            end    = int(metadata["SubjectEnd"])

            actual_species = species_from_id(seq_id)
            matches = [f for f in os.listdir(args.genome_dir)
                       if f.startswith(actual_species) and f.endswith(".fasta")]
            if not matches:
                print(f"Warning: No genome file starting with '{actual_species}' in {args.genome_dir}. Skipping.")
                skipped += 1
                continue
            genome_filepath = os.path.join(args.genome_dir, matches[0])

            upstream_seq, gene_seq, downstream_seq = extract_region(
                genome_filepath, contig, start, end, upstream_bp, downstream_bp
            )

            if upstream_seq is None:
                skipped += 1
                continue

            if merge:
                full_seq = upstream_seq + gene_seq + downstream_seq
                write_record(merged_fh, header, seq_id, dict(metadata), "merged", full_seq)
            else:
                if up_fh is not None:
                    seq = (upstream_seq + gene_seq) if not flank_only else upstream_seq
                    if seq:
                        write_record(up_fh, header, seq_id, dict(metadata), "upstream", seq)

                if dn_fh is not None:
                    seq = (gene_seq + downstream_seq) if not flank_only else downstream_seq
                    if seq:
                        write_record(dn_fh, header, seq_id, dict(metadata), "downstream", seq)

    for fh in (up_fh, dn_fh, merged_fh):
        if fh is not None:
            fh.close()

    if skipped:
        print(f"Skipped {skipped} record(s) due to missing data.")

    written = []
    if merge:
        written.append("merged.fasta")
    else:
        if upstream_bp   > 0: written.append("upstream.fasta")
        if downstream_bp > 0: written.append("downstream.fasta")
    print(f"Neighbour sequences written to: {', '.join(written)}")


if __name__ == "__main__":
    main()
