#!/usr/bin/env python
"""
Re-wrap a FASTA file (default 60 nt/aa per line).

Usage
-----
python format_fasta.py input.fasta output.fasta
"""

import os
import sys
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter


def reformat_fasta(input_file: str, output_file: str, line_length: int = 60) -> None:
    # Make sure the destination directory exists
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)

    # Read records, write with specified wrap length
    with open(input_file) as in_handle, open(output_file, "w") as out_handle:
        writer = FastaWriter(out_handle, wrap=line_length)
        writer.write_file(SeqIO.parse(in_handle, "fasta"))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python format_fasta.py input.fasta output.fasta")

    in_fasta, out_fasta = sys.argv[1], sys.argv[2]
    reformat_fasta(in_fasta, out_fasta)
    print(f"✔ Reformatted FASTA saved to {out_fasta}")
