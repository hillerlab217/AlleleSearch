# Allele Search Pipeline

A Nextflow pipeline for identifying protein allele variants from genomic sequences. It searches for query protein/DNA sequences across multiple genomes using BLAST, identifies open reading frames, translates them to protein, and clusters identical proteins into allele groups with variant annotation.

## Pipeline Overview

```
Query sequences + Genomes
        │
        ▼
   BLAST Search ──► Parse & Filter hits
        │
        ▼
   ORF Finding ──► Check start/stop codons ──► Extend incomplete sequences
        │
        ▼
   Translation ──► All 3 reading frames, keep longest ORF
        │
        ▼
   Allele Analysis ──► Group identical proteins, identify variants
        │
        ▼
   Alignment & Visualization ──► MAFFT, WebLogo, T-Coffee
```

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
6. **Quality & Visualization** — Filters alleles by species representation, generates MAFFT alignments, WebLogo sequence logos, and optionally T-Coffee alignments.

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
├── allele_info/
│   └── <query_name>/
│       ├── query_reference/    # Alleles compared to original query
│       └── self_reference/     # Alleles identified by self-comparison
├── allele_seq_filtered/        # Alleles passing species threshold
├── allele_seq/                 # MAFFT-aligned allele sequences
├── logo/                       # WebLogo sequence logos (.eps)
├── allele_tcoffee/             # T-Coffee alignments (if enabled)
├── timeline.html               # Nextflow execution timeline
└── report.html                 # Nextflow execution report
```

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
    └── protein_allele_variant.py     # Allele grouping and variant caller
```
