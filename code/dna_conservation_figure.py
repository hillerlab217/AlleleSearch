#!/usr/bin/env python3
"""
dna_conservation_figure.py

Given a MAFFT-aligned DNA FASTA, compute per-column nucleotide conservation
scores and draw a line figure (Conservation score vs. Genomic position).

Conservation score at each column = frequency of the most common nucleotide
among non-gap characters. Columns where the consensus is a gap are excluded;
remaining positions are numbered sequentially (1 to N).

Sequences are deduplicated by genome (the genome ID extracted from the FASTA
header after "vs_") before scoring, so multiple overlapping BLAST hits from
the same genome cannot inflate conservation.

Usage
-----
python dna_conservation_figure.py ALIGNED.fasta OUTPUT QUERY_NAME [--format pdf|png|svg]
"""

import sys
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO


# ── helpers ──────────────────────────────────────────────────────────────────

def genome_id_from_header(header: str) -> str:
    """Extract genome identifier from a FASTA header ID.

    Headers follow the pipeline convention:
        >query_vs_genome_name|contig=...|...

    Returns the part after 'vs_' (the genome name).  If 'vs_' is not present,
    the full ID is returned unchanged so deduplication still works.
    """
    seq_id = header.split('|')[0]
    if 'vs_' in seq_id:
        return seq_id.split('vs_', 1)[1]
    return seq_id


def deduplicate_by_genome(records):
    """Keep only the first BLAST hit per genome.

    Multiple overlapping BLAST hits from the same genome would bias
    conservation scores toward 1.0 (identical sequences counted multiple times).
    """
    seen = set()
    kept = []
    for rec in records:
        gid = genome_id_from_header(rec.id)
        if gid not in seen:
            seen.add(gid)
            kept.append(rec)
    return kept


def compute_conservation(records):
    """Return a per-column conservation score array.

    Score = count of most common non-gap nucleotide / total non-gap count.
    All-gap columns receive score 0.0 (no division attempted).
    """
    n_cols = len(records[0].seq)
    scores = np.zeros(n_cols)
    for i in range(n_cols):
        col = [str(rec.seq[i]).upper() for rec in records]
        non_gap = [c for c in col if c not in ('-', 'N', 'n')]
        if not non_gap:          # Fix #2: guard against division by zero
            scores[i] = 0.0
            continue
        counts = Counter(non_gap)
        scores[i] = counts.most_common(1)[0][1] / len(non_gap)
    return scores


def compute_consensus(records):
    """Plurality-rule consensus; gaps treated as absent."""
    n_cols = len(records[0].seq)
    consensus = []
    for i in range(n_cols):
        col = [str(rec.seq[i]).upper() for rec in records]
        non_gap = [c for c in col if c != '-']
        if not non_gap:
            consensus.append('-')
        else:
            consensus.append(Counter(non_gap).most_common(1)[0][0])
    return ''.join(consensus)


# ── drawing ───────────────────────────────────────────────────────────────────

def draw_conservation_figure(positions, scores, query_name, output_path):
    fig, ax = plt.subplots(figsize=(12, 4))

    ax.plot(positions, scores, color='#4f81bd', linewidth=0.9, label=query_name)

    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(positions[0], positions[-1])
    ax.set_xlabel("Genomic position", fontsize=11)
    ax.set_ylabel("Conservation score", fontsize=11)
    ax.legend(frameon=False, fontsize=9)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    suffix = Path(output_path).suffix.lower()
    dpi = 300 if suffix == '.png' else 150
    fig.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Wrote {output_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Draw a nucleotide conservation score figure from a MAFFT-aligned FASTA."
    )
    parser.add_argument("aligned_fasta", help="MAFFT-aligned DNA FASTA file")
    parser.add_argument("output",        help="Output figure path")
    parser.add_argument("query_name",    help="Query name (used as plot label)")
    parser.add_argument("--format", default="pdf",
                        choices=["pdf", "png", "svg", "eps", "jpeg"],
                        help="Output format (default: pdf)")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.aligned_fasta, "fasta"))

    if len(records) < 2:
        sys.exit(
            f"ERROR: Need at least 2 sequences to compute conservation; "
            f"found {len(records)} in {args.aligned_fasta}"
        )

    # Fix #3: deduplicate — one sequence per genome
    records = deduplicate_by_genome(records)
    if len(records) < 2:
        sys.exit(
            f"ERROR: After deduplication by genome, fewer than 2 unique genomes remain. "
            f"Cannot compute conservation."
        )

    scores    = compute_conservation(records)
    consensus = compute_consensus(records)

    # Keep only positions where consensus is not a gap; number them 1..N
    positions       = []
    filtered_scores = []
    pos = 0
    for i, c in enumerate(consensus):
        if c != '-':
            pos += 1
            positions.append(pos)
            filtered_scores.append(scores[i])

    if not positions:
        sys.exit("ERROR: All alignment columns are gap-only — nothing to plot.")

    draw_conservation_figure(positions, filtered_scores, args.query_name, args.output)


if __name__ == "__main__":
    main()
