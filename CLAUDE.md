# Allele Search Pipeline — CLAUDE.md

## Project Context

This is a bioinformatics project using Python. When debugging, always check input data constraints (e.g., sequence lengths, file formats) before diving into code logic.

## Workflow

When working through a checklist or TODO list, always re-read the checklist file after completing each major item to track progress and avoid missing steps.

---

## Project Purpose

Nextflow DSL2 pipeline for discovering protein allele variants across large collections of bacterial/microbial genomes. Given one or more query DNA/protein sequences and a directory of genome FASTA files, it:

1. BLASTs each query against all genomes
2. Validates and extends open reading frames (start/stop codon checking)
3. Translates DNA hits to proteins
4. Clusters identical proteins into allele groups
5. Calls amino acid variants relative to a reference
6. Computes a consensus reference sequence from all alleles
7. Produces colored alignment figures and WebLogo sequence logos

**Domain context:** "Allele" here means a unique protein variant found in ≥1 sequenced strain. Variant positions are formatted as `90:K→T` (position : reference AA → observed AA). The most frequent allele is named `allele_1`.

---

## Tech Stack

| Tool | Version | Role |
|---|---|---|
| Nextflow | DSL2 | Workflow orchestration |
| Python | 3.9 | All custom scripts (`code/`) |
| BioPython | 1.85 | Sequence I/O and alignment |
| BLAST+ | 2.16.0 | `blastn` / `tblastn` searches |
| MAFFT | latest | Multiple sequence alignment |
| WebLogo | latest | Sequence logo generation |
| Matplotlib | 3.9.4 | Alignment figure rendering |

**Conda environment:** `allele_search` (defined in `environment.yml`)

---

## Environment Setup

```bash
# One-time setup
conda env create -f environment.yml
conda activate allele_search

# Verify pipeline
nextflow run main.nf --help
```

`nextflow.config` hard-codes the conda env path in `beforeScript`:
```
export PATH=/home/xinruoz/miniconda3/envs/allele_search/bin:$PATH
```
**Update this path in `nextflow.config` if the conda prefix changes.**

---

## Running the Pipeline

### Full pipeline (most common)
```bash
nextflow run main.nf \
  --query_dir ./data \
  --genome_dir /path/to/genomes \
  --outdir ./results \
  --sequence_identity 98 \
  --length_percentage 98 \
  -resume
```

### Step-specific modes
Each mode reads from and writes to `--outdir`. Prior steps must have completed first.

```bash
# Re-run BLAST only (creates/rebuilds blast DBs, runs search, parses, filters)
nextflow run main.nf --blast --query_dir ./data --genome_dir /genomes --outdir ./results

# Re-filter with different thresholds (reuses existing blast_summary/ CSVs)
nextflow run main.nf --filter --query_dir ./data --outdir ./results --sequence_identity 95

# Re-run ORF finding + translation (reuses existing blast_summary_filtered/ CSVs)
nextflow run main.nf --orf --genome_dir /genomes --outdir ./results

# Re-run visualization only (reuses existing allele_info/ and allele_seq_filtered/)
nextflow run main.nf --visual --outdir ./results --weblogo_format png

# DNA conservation figure — run own BLAST (blastn) from scratch
nextflow run main.nf --dna_variation --query_dir ./data --genome_dir /genomes --outdir ./results

# DNA conservation figure — reuse existing dna_variation/blast_summary_filtered/ CSVs
nextflow run main.nf --dna_variation --outdir ./results -resume
```

Alternative syntax: `--step blast|filter|orf|visual|dna_variation`

### Test run
```bash
nextflow run main.nf \
  --query_dir ./test_query \
  --genome_dir /path/to/test_genomes \
  --outdir ./results/test
```
Test query: `test_query/murM.fa`

---

## Directory & Module Map

```
main.nf                          # Workflow entrypoint; routes step modes
nextflow.config                  # All parameters with defaults
environment.yml                  # Conda dependencies
modules/
  blast.nf                       # COLLECT_ALL_BLAST_DBS, BLAST_SEARCH, PARSE_BLAST
  filter.nf                      # FILTER_BLAST
  orf.nf                         # EXTRACT_DNA, CSV_TO_FASTA, CODON_CHECK,
                                 #   EXTEND_SEQUENCES, COLLECT_STANDARD
  translate.nf                   # TRANSLATE (DNA → protein, longest ORF)
  allele.nf                      # ALLELE_SELF_REF, ALLELE_CONSENSUS_REF
  quality.nf                     # MAFFT_FOR_CONSENSUS, COMPUTE_CONSENSUS,
                                 #   FILTER_ALLELES, MAFFT_ALIGN, WEBLOGO,
                                 #   ALIGNMENT_FIGURE, ALIGNMENT_FIGURE_REF
  dna_variation.nf               # BLAST_SEARCH_DNA, PARSE_BLAST_DNA, FILTER_BLAST_DNA,
                                 #   MAFFT_DNA, DNA_CONSERVATION_FIGURE
code/
  header_utils.py                # Canonical FASTA header parser/builder (shared utility)
  blastParser.py                 # Converts raw BLAST .out → CSV
  filter_blast.py                # Filters parsed BLAST CSV by identity/length/e-value
  csv_to_fasta.py                # CSV → FASTA with metadata headers
  start_stop_codon_dna_seq.py   # Classifies sequences by codon validity; trims
  up_down_stream.py              # Extends incomplete sequences from genome files
  extract_dna_from_blast_csv.py  # Replaces tblastn protein seqs with DNA coords
  translate.py                   # DNA → protein in 3 frames, returns longest ORF
  compute_consensus.py           # Computes consensus sequence from MAFFT alignment
  protein_allele_variant.py      # Core: groups identical proteins, calls variants (BLOSUM62)
  alignment_figure.py            # Renders colored AA alignment PDFs
  dna_conservation_figure.py     # Renders per-position DNA conservation score line plots (--dna_variation)
  format_fasta.py                # Utility: FASTA formatting helpers
  highlight_variant_pos_protein_structure.py  # Auxiliary: highlights variant positions on protein structures (not called by pipeline)
  pyMol.py                       # Auxiliary: PyMOL visualization helper (not called by pipeline)
data/                            # Example query FASTA files (insertion regions)
test_query/                      # Test data (murM.fa)
results/                         # Pipeline outputs (gitignored)
```

---

## Full Pipeline Step-by-Step

| Step | Process | Description |
|---|---|---|
| 1 | `COLLECT_ALL_BLAST_DBS` | Auto-build missing BLAST nucleotide DBs for all genome FASAs |
| 2 | `BLAST_SEARCH` | Run blastn/tblastn in batches of 100 genomes per task |
| 3 | `PARSE_BLAST` | Merge batch `.out` files → single CSV per query |
| 4 | `FILTER_BLAST` | Apply identity/length/e-value cutoffs |
| 5 | `EXTRACT_DNA` | tblastn only: replace protein hits with DNA coordinates from genome |
| 6 | `CSV_TO_FASTA` | Convert filtered CSV → FASTA with metadata headers |
| 7 | `CODON_CHECK_10` | Classify sequences by start/stop codons (10-codon window) |
| 8 | `EXTEND_SEQUENCES` | Fetch upstream/downstream sequence from genome for incomplete ORFs |
| 9 | `CODON_CHECK_20` | Re-classify extended sequences (20-codon window) |
| 10 | `COLLECT_STANDARD` | Merge sequences that passed both codon rounds |
| 11 | `TRANSLATE` | DNA → protein (3-frame, longest ORF wins) |
| 12 | `ALLELE_SELF_REF` | Group identical proteins into alleles; call variants vs. most frequent allele |
| 13 | `MAFFT_FOR_CONSENSUS` | Align all allele sequences |
| 14 | `COMPUTE_CONSENSUS` | Derive consensus reference from the alignment |
| 15 | `ALLELE_CONSENSUS_REF` | Re-call variants vs. consensus reference |
| 16 | `FILTER_ALLELES` | Remove alleles below `--species_threshold` % |
| 17 | `ALIGNMENT_FIGURE` | PDF alignment figure for self-reference (only when `--self_reference_figures`) |
| 18 | `ALIGNMENT_FIGURE_REF` | PDF alignment figure for consensus-reference (always generated; consensus row labeled `{query_name}_consensus`) |
| 19 | `MAFFT_ALIGN` | Align filtered alleles for WebLogo |
| 20 | `WEBLOGO` | Generate sequence logo (weighted by strain frequency) |

---

## Key Conventions & Design Decisions

### FASTA header format (internal)
All scripts communicate metadata through structured FASTA headers, parsed and built exclusively via `header_utils.py`. **Do not change this format without updating `header_utils.py`.**
```
>ID|contig=NAME|E-value=VALUE|AlignLen=LEN|SeqIdentity=%|SubjectStart=POS|SubjectEnd=POS
```

### Strand detection
**SubjectStart > SubjectEnd signals reverse strand.** No explicit strand column exists. Used in `up_down_stream.py`, `extract_dna_from_blast_csv.py`, and `start_stop_codon_dna_seq.py`.

### Species name extraction
Parsed from BLAST sequence IDs by splitting on `"vs_"` (implemented in `header_utils.species_from_id()`). BLAST output files follow `{query}_vs_{genome_name}`. Changing this naming convention breaks allele species counting.

### Two allele reference modes
Both allele info CSVs are always generated. Alignment figures differ:
- `self_reference/`: variants called against the most frequent allele. Alignment figures only generated when `--self_reference_figures` is set (default: off).
- `consensus_reference/`: variants called against a computed consensus of all alleles. Alignment figures always generated; the consensus reference row is labeled `{query_name}_consensus` (e.g. `murM_consensus`).

### BLAST batching
Genomes are batched 100 per task (configurable: `--blast_batch_size`). Each batch produces a separate `.out` file; `PARSE_BLAST` collects all batches per query with `groupTuple()`.

### BLAST DB creation
`COLLECT_ALL_BLAST_DBS` always auto-checks and builds any missing `.nin`-indexed databases before BLAST runs. No manual flag needed.

### Codon check logic
- Round 1 (`CODON_CHECK_10`): 10-codon window; searches all 3 reading frames
- Round 2 (`CODON_CHECK_20`): 20-codon window on genome-extended sequences
- Valid start codons: `ATG`, `GTG`, `TTG`
- Valid stop codons: `TAA`, `TAG`, `TGA`
- Sequences with internal in-frame stop codons (`more_than_one_stop`) are discarded
- Classified outputs: `standard`, `no_stop`, `no_start`, `no_stop_no_start`, `stop_not_at_end`, `no_start_stop_not_at_end`, `more_than_one_stop`

### Consensus sequence
Computed by MAFFT-aligning all alleles, then `compute_consensus.py` picks the plurality AA at each position. Used as the reference for `ALLELE_CONSENSUS_REF` and `ALIGNMENT_FIGURE_REF`.

### WebLogo weighting
Before WebLogo runs, sequences are duplicated to weight the logo by strain frequency. The awk in `quality.nf` first uses `species_number=N` (raw strain count) if present; otherwise it rounds `species_percentage` to the nearest int. Headers must contain at least one of these fields.

### Amino acid color scheme (alignment figures)
Physicochemical groups in `alignment_figure.py`:
- Orange `#f4b183`: Hydrophobic (A G I L M V P)
- Blue `#4f81bd`: Positive (R K H)
- Green `#9bbb59`: Polar (N C Q S T)
- Red `#e26b6b`: Negative (D E)
- Yellow `#fff475`: Aromatic (F W Y)

---

## Parameters Reference

| Parameter | Default | Description |
|---|---|---|
| `--query_dir` | required | Directory with query `.fasta/.fa/.fna/.ffn` files |
| `--genome_dir` | required | Directory with genome FASTA files |
| `--outdir` | `./results` | Output directory |
| `--blast_type` | `blastn` | `blastn` or `tblastn` |
| `--blast_batch_size` | `100` | Genomes per BLAST task (not in `--help`, set in `nextflow.config`) |
| `--max_evalue` | `0.00001` | E-value cutoff |
| `--sequence_identity` | `80` | Minimum % sequence identity |
| `--length_percentage` | `80` | Minimum alignment length as % of query length |
| `--extension_codons` | `10` | Codons to extend incomplete ORFs |
| `--species_threshold` | `1` | Minimum species % for allele filtering / figures |
| `--self_reference_figures` | `false` | Also generate self-reference alignment figures (default: only consensus-reference figures are generated) |
| `--dna_variation_figure_format` | `pdf` | Output format for the DNA conservation figure: `pdf`, `png`, `svg`, `eps`, `jpeg` |
| `--weblogo_format` | `eps` | `eps`, `png`, `pdf`, `svg`, `jpeg` |
| `--weblogo_units` | `probability` | `probability`, `bits`, `nats` |
| `--weblogo_color_scheme` | `chemistry` | `chemistry`, `hydrophobicity`, `charge`, `classic` |
| `--weblogo_size` | `medium` | `small`, `medium`, `large` |
| `--weblogo_resolution` | `96` | DPI |

---

## Common Gotchas

- **`--filter` mode requires `--query_dir`** — needed to calculate minimum alignment length from query length.
- **`--filter` mode requires prior `--blast` run** — needs `{outdir}/blast_summary/*.csv` to exist.
- **`--orf` mode requires `--genome_dir`** — needed for upstream/downstream sequence extension.
- **`--visual` mode requires prior full run** — needs `{outdir}/allele_seq_filtered/` and `{outdir}/allele_info/`.
- **`--blast_type tblastn` triggers `EXTRACT_DNA`** — swaps protein subject sequences for the underlying DNA coordinates before ORF checking.
- **MAFFT is skipped silently** if fewer than 2 sequences remain after `FILTER_ALLELES`. Search for "Skipping MAFFT" in Nextflow logs (the message does not use the word WARNING).
- **`COLLECT_STANDARD` emits an optional output** — if no sequences survive both codon rounds, downstream steps are skipped automatically (no crash).
- **`-resume` reuses cached tasks** — omit it to force a full re-run. Step modes always overwrite their output subdirectories regardless of `-resume`.
- **Nextflow work directory grows large** — run `nextflow clean -f` periodically.

---

## Development Workflow

### Adding a new parameter
1. Add to `params {}` in `nextflow.config` with a sensible default
2. Reference as `${params.param_name}` in the relevant module `script:` block
3. Add to the `--help` message in `main.nf`
4. Document in `README.md`

### Adding a new pipeline step
1. Create or edit the relevant `modules/*.nf` file
2. Import the process at the top of `main.nf`
3. Wire it into the appropriate workflow block(s) (`full` and any relevant step modes)
4. Update channel routing if downstream steps depend on its output

### Modifying Python scripts
Scripts are called from Nextflow `script:` blocks via:
```
python ${params.code_dir}/script_name.py [args]
```
`code_dir` defaults to `${projectDir}/code`. All scripts use `argparse`. Header parsing/building must go through `header_utils.py`.

### Testing a step change
```bash
# Re-run only the affected step; -resume reuses upstream cached work
nextflow run main.nf --visual --outdir ./results/test -resume
```

---

## Output Interpretation

### Key outputs: `results/allele_info/<query_name>/`

| Path | Description |
|---|---|
| `self_reference/allele_info_no_strain_name.csv` | Alleles ranked by frequency; variants vs. most frequent allele |
| `self_reference/allele_info_with_strain_name.csv` | Same, with strain names listed per allele |
| `self_reference/allele_sequence_output.fasta` | Protein FASTA, one entry per allele |
| `consensus_reference/allele_info_no_strain_name.csv` | Alleles ranked by frequency; variants vs. consensus |
| `consensus_reference/allele_sequence_output.fasta` | Protein FASTA, one entry per allele (same format, consensus-reference run) |
| `consensus_reference/<query>_consensus.fasta` | The computed consensus reference sequence |

CSV columns: `allele_name`, `species_count`, `species_percentage`, `variant_position_count`, `variant_percentage`, `variant_positions(position:reference→changed)`, `allele_sequence`

### Alignment figures: `results/alignment_figures/<query>/{self_reference,consensus_reference}/`
- `*_alignment.pdf` — only variant positions shown
- `*_alignment_gap.pdf` — all positions between first and last variant

### Logos: `results/logo/<query_name>_logo.<format>`

### Intermediate outputs
| Directory | Contents |
|---|---|
| `blast_raw/` | Raw BLAST `.out` files (one per query-batch) |
| `blast_summary/` | Parsed CSVs (one per query, all hits) |
| `blast_summary_filtered/` | Filtered CSVs after identity/length/e-value cutoffs |
| `standard_dna_seq/` | DNA sequences passing both codon-check rounds |
| `allele_seq_filtered/` | Per-allele protein FASTAs above species threshold |
| `allele_seq/` | MAFFT-aligned filtered allele FASTAs |

---

## Related Documentation

- `README.md` — User-facing installation and usage guide
- `run_example.txt` — Real command examples from previous analysis runs
- `task.md` — Development task history (all phases complete)
- `flowchart.mermaid` — Visual overview of the pipeline DAG
