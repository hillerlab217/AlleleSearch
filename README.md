# Allele Search Pipeline

A Nextflow pipeline for identifying protein allele variants from genomic sequences. It searches for query protein/DNA sequences across multiple genomes using BLAST, screens proteins according to the quality of BLAST results, and ensures start and stop codons are present, no stop codons are in between, identifies open reading frames, translates them into proteins, and clusters identical proteins into allele groups with variant annotations.

## Pipeline Overview

```
Query sequences + Genomes
        │
        ▼
   BLAST Search ──► Blast, Parse Filter hits
        │
        ▼
   ORF Finding ──► Check start/stop codons ──► Extend incomplete sequences
        │
        ▼
   Translation ──► All 3 reading frames, keep longest ORF
        │
        ▼
   Allele Alignment ──► Group identical proteins, identify variants
        │
        ▼
   Visualization ──► Alignment Figures, WebLogo, T-Coffee
```

## Core outputs



## Installation

The pipeline requires macOS or Linux terminal to run, miniconda/conda to help install and manage software packages.

### Step 1: Download This Pipeline

**Option A: Download as ZIP**
1. Go to the repository page on GitHub
2. Click the green "Code" button
3. Click "Download ZIP"
4. Extract the ZIP file and navigate to the folder in terminal

**Option B: Clone with Git**
```bash
git clone <repository-url> allele_search
cd allele_search
```

### Step 2: Set Up the Environment

```bash
conda env create -f environment.yml
conda activate allele_search
```

**Note:** Run `conda activate allele_search` every time you open a new terminal before using the pipeline.

### Step 3: Verify Installation

```bash
nextflow run main.nf --help
```

You should see a help message listing all available options.

### Steps in Detail

1. **BLAST Search** — Runs each query against each genome (blastn or tblastn). Parses results to CSV.
2. **Filtering** — Keeps hits passing e-value, similarity %, and alignment length thresholds.
3. **ORF Finding** — Checks for start/stop codons in a 10-codon window. Sequences missing codons are extended using the reference genome, then re-checked with a 20-codon window.
4. **Translation** — Translates DNA to protein in all 3 frames. The longest ORF (M to stop) is kept.
5. **Allele Identification** — Groups identical protein sequences into alleles. Variant positions are identified by aligning each allele to a reference using BLOSUM62.
6. **Quality & Visualization** — Filters alleles by species representation, generates coloured variant alignment figures, MAFFT alignments, WebLogo sequence logos, and optionally T-Coffee alignments.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Conda (environments are managed automatically via `environment.yml`)
- BLAST+
- Python 3 with BioPython
- MAFFT
- WebLogo
- T-Coffee (optional)

## Usage

```bash
nextflow run main.nf \
  --query_dir <path/to/query/fastas> \
  --genome_dir <path/to/genome/fastas>
```

Example with all optional parameters:

```bash
nextflow run main.nf \
  --query_dir ./queries \
  --genome_dir ./genomes \
  --outdir ./my_results \
  --blast_type tblastn \
  --make_blastdb true \
  --max_evalue 1e-10 \
  --similarity 90 \
  --length_percentage 85 \
  --extension_codons 20 \
  --species_threshold 5 \
  --skip_tcoffee true
```

Show all options:

```bash
nextflow run main.nf --help
```

### Required Parameters

| Parameter | Description |
|---|---|
| `--query_dir` | Directory containing query `.fasta` files |
| `--genome_dir` | Directory containing genome files (`.fasta`, `.fa`, `.fna`, `.ffn`). Also used for sequence extension |

### Optional Parameters

| Parameter | Default | Description |
|---|---|---|
| `--outdir` | `./results` | Output directory |
| `--blast_type` | `blastn` | BLAST program (`blastn` or `tblastn`) |
| `--make_blastdb` | `false` | Build new BLAST databases. If `false`, uses existing `.n*` DB files |
| `--max_evalue` | `0.00001` | E-value threshold |
| `--similarity` | `80` | Minimum similarity percentage |
| `--length_percentage` | `80` | Minimum alignment length as percentage of query length |
| `--extension_codons` | `10` | Number of codons to extend incomplete sequences |
| `--species_threshold` | `1` | Minimum species percentage to keep an allele |
| `--skip_tcoffee` | `false` | Skip T-Coffee alignment step |



## Output

Results are written to `--outdir` (default `./results/`):

```
results/
├── blast_summary/              # Raw BLAST results (CSV)
├── blast_summary_filtered/     # Filtered BLAST results (CSV)
├── standard_dna_seq/           # DNA sequences with valid start/stop codons
├── protein_seq/
│   ├── inframe/                # Longest ORF protein sequences
│   └── others/                 # Alternative reading frames
├── allele_info/                # KEY OUTPUT: per-query allele variant analysis
│   └── <query_name>/
│       ├── query_reference/    # Variants compared to the original query sequence
│       │   ├── allele_info_with_strain_name.csv
│       │   ├── allele_info_no_strain_name.csv
│       │   └── allele_sequence_output.fasta
│       └── self_reference/     # Variants compared to the most abundant allele found
│           ├── allele_info_with_strain_name.csv
│           ├── allele_info_no_strain_name.csv
│           └── allele_sequence_output.fasta
├── alignment_figures/           # Coloured variant alignment figures (PDF)
│   └── <query_name>/
│       ├── self_reference/     # Figures using most abundant allele as reference
│       └── query_reference/    # Figures using original query as reference
├── allele_seq_filtered/        # Alleles passing species threshold
├── allele_seq/                 # MAFFT-aligned allele sequences
├── logo/                       # WebLogo sequence logos (.eps)
├── allele_tcoffee/             # T-Coffee alignments (if enabled)
├── timeline.html               # Nextflow execution timeline
└── report.html                 # Nextflow execution report
```

### `allele_info/` — Key Output Details

Each CSV row represents one allele (a group of identical protein sequences) and includes:

| Column | Description |
|--------|-------------|
| `allele_name` | Allele ID, numbered by frequency (1 = most common) |
| `species_count` | Number of species carrying this allele |
| `species_percentage` | Percentage of all searched species with this allele |
| `variant_position_count` | Number of amino acid positions differing from reference |
| `variant_percentage` | Percentage of positions that are variants |
| `variant_positions` | Specific amino acid changes, e.g. `90:K→T; 93:Q→K` |
| `allele_sequence` | Full protein sequence of the allele |
| `strain_names` | Strain list (only in `_with_strain_name.csv`) |

- **`query_reference/`** uses the original query as the reference for variant calling.
- **`self_reference/`** uses the most abundant allele found as the reference.
- **`allele_sequence_output.fasta`** contains the protein sequence of each allele, with species count in the header.

## Project Structure

```
├── main.nf              # Main pipeline workflow
├── nextflow.config       # Pipeline configuration
├── modules/
│   ├── blast.nf          # BLAST database creation and search
│   ├── filter.nf         # Hit filtering by quality metrics
│   ├── orf.nf            # ORF finding, codon checking, sequence extension
│   ├── translate.nf      # DNA to protein translation
│   ├── allele.nf         # Allele grouping and variant identification
│   └── quality.nf        # Alignment, filtering, and visualization
└── code/
    ├── blastParser.py                # BLAST output parser
    ├── csv_to_fasta.py               # CSV to FASTA converter
    ├── start_stop_codon_dna_seq.py   # Start/stop codon classifier
    ├── up_down_stream.py             # Sequence extension from reference genomes
    ├── translate.py                  # DNA to protein translator
    ├── protein_allele_variant.py     # Allele grouping and variant caller
    └── alignment_figure.py           # Variant alignment figure generator
```

