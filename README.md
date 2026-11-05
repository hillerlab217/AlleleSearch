# Allele Search Pipeline

A Nextflow pipeline for identifying protein allele variants from genomic sequences. It searches for query protein/DNA sequences across multiple genomes using BLAST, screens proteins according to the quality of BLAST results, and ensures start and stop codons are present, no stop codons are in between, identifies open reading frames, translates them into proteins, and clusters identical proteins into allele groups with variant annotations.

## Pipeline Overview

```
Query sequences + Genomes
        │
        ▼
   BLAST Search ──► Blast, Parse, Filter hits
        │
        ├── [--mode protein (default)] ──────────────────────────────┐
        │    ORF Finding ──► Check start/stop codons                 │
        │         │                                                  │
        │         ▼                                                  │
        │    Translation ──► Frame 0 protein, M-start enforced       │
        │         │                                                  │
        └── [--mode dna] ──► (ORF + translation skipped) ──────────►─┘
                                                                     │
        ▼                                                            │
   Allele Analysis ──► Group identical sequences, call variants ◄───┘
        │
        ▼
   Visualization ──► Alignment Figures, WebLogo
```

## Core outputs

| Path | Description |
|---|---|
| `allele_info/<query>/consensus_reference/allele_info_no_strain_name.csv` | Alleles ranked by frequency; variants called vs. a majority-vote consensus sequence |
| `allele_info/<query>/self_reference/allele_info_no_strain_name.csv` | Same, with variants called vs. the most abundant allele |
| `alignment_figures/<query>/consensus_reference/` | Coloured variant alignment PDFs (consensus reference) |
| `logo/<query>_logo.<format>` | WebLogo sequence logo weighted by strain frequency |

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

After creating the environment, update the conda path in `nextflow.config` to match your local miniconda prefix:
```
process.beforeScript = 'export PATH=/path/to/miniconda3/envs/allele_search/bin:$PATH'
```

### Step 3: Verify Installation

```bash
nextflow run main.nf --help
```

You should see a help message listing all available options.

### Steps in Detail

1. **BLAST Search** — Runs each query against each genome (blastn or tblastn). Parses results to CSV.
2. **Filtering** — Keeps hits passing e-value, sequence identity %, and alignment length thresholds.
3. **ORF Finding** *(protein mode only)* — Checks for start/stop codons in a 10-codon window. Sequences missing codons are extended using the reference genome, then re-checked with a 20-codon window. Extended sequences shorter than 5 amino acids are discarded.
4. **Translation** *(protein mode only)* — Translates the validated reading frame (Frame 0) to protein. The start codon is replaced with ATG, so all output proteins begin with Met regardless of the original start codon (ATG, GTG, or TTG).
5. **Allele Identification** — Groups identical sequences into alleles. In protein mode, variant positions are identified using BLOSUM62 alignment. In DNA mode, nucleotide-level variants are called (e.g. `90:A→T`).
6. **Quality & Visualization** — Filters alleles by species representation, generates coloured variant alignment figures, MAFFT alignments, and WebLogo sequence logos.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2) — install separately: `curl -s https://get.nextflow.io | bash`
- Conda (environments are managed automatically via `environment.yml`)
- BLAST+
- Python 3 with BioPython
- MAFFT
- WebLogo

> All dependencies except Nextflow are installed automatically by the conda environment.

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
  --sequence_identity 90 \
  --length_percentage 85 \
  --extension_codons 20 \
  --species_threshold 5
```

ATG-only starts (stricter), exclude TAG as stop codon:

```bash
# ATG-only starts (stricter), exclude TAG as stop codon
nextflow run main.nf --query_dir ./data --genome_dir /genomes --outdir ./results \
  --start_codons ATG \
  --stop_codons TAA,TGA
```

> **Note:** When using non-ATG start codons (e.g. `GTG`, `TTG`), the start codon is
> translated as Met (M) in the output protein, which reflects the biological reality of
> fMet initiation at all start codons.

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
| `--mode` | `protein` | Sequence mode: `protein` (default) or `dna`. In `dna` mode, ORF finding and translation are skipped; BLAST hits are grouped and annotated directly as nucleotide alleles |
| `--blast_type` | `blastn` | BLAST program (`blastn` or `tblastn`) |
| `--make_blastdb` | `false` | Build new BLAST databases. If `false`, uses existing `.n*` DB files |
| `--max_evalue` | `0.00001` | E-value threshold |
| `--sequence_identity` | `80` | Minimum sequence identity percentage |
| `--length_percentage` | `80` | Minimum alignment length as percentage of query length |
| `--extension_codons` | `10` | Number of codons to extend incomplete sequences |
| `--start_codons` | `ATG,GTG,TTG` | Comma-separated valid start codons |
| `--stop_codons` | `TAA,TGA,TAG` | Comma-separated valid stop codons |
| `--species_threshold` | `1` | Minimum species percentage to keep an allele |
| `--self_reference_figures` | `false` | Also generate self-reference alignment figures; by default only consensus-reference figures are generated |
| `--blast_batch_size` | `100` | Number of genomes processed per BLAST task |

### Neighbour Extraction Parameters

Extract flanking sequence (upstream/downstream of the gene) from the reference genomes. These can be used standalone with `--neighbour` or added to a full run.

| Parameter | Default | Description |
|---|---|---|
| `--upstream` | `0` | bp upstream of gene to extract (0 = skip) |
| `--downstream` | `0` | bp downstream of gene to extract (0 = skip) |
| `--flank_only` | `false` | Exclude the gene itself; output flanking regions only |
| `--neighbour_merge` | `false` | Merge upstream + gene + downstream into one FASTA (requires gene to be included) |

### BLAST Tuning Parameters

These parameters are passed directly to the BLAST command. All default to unset (BLAST uses its own defaults).

**Shared — apply to both `blastn` and `tblastn`:**

| Parameter | Description |
|---|---|
| `--blast_word_size` | Word size for seeding (`-word_size`). Min 4 for blastn, min 2 for tblastn |
| `--blast_gapopen` | Gap open penalty (`-gapopen`) |
| `--blast_gapextend` | Gap extend penalty (`-gapextend`) |

**tblastn only — ignored when `--blast_type blastn`:**

| Parameter | Description |
|---|---|
| `--tblastn_matrix` | Substitution matrix (`-matrix`), e.g. `PAM30`, `BLOSUM45` (default: `BLOSUM62`) |
| `--tblastn_threshold` | Minimum word score threshold (`-threshold`) |
| `--tblastn_comp_based` | Composition-based statistics (`-comp_based_stats`): `0`=off, `1`, `2` (default) |
| `--tblastn_seg` | SEG low-complexity filter (`-seg`): `yes` or `no` |

**Example — tblastn tuned for short/divergent peptide queries:**

```bash
nextflow run main.nf \
  --query_dir ./queries \
  --genome_dir ./genomes \
  --outdir ./results/my_analysis \
  --blast_type tblastn \
  --max_evalue 200000 \
  --tblastn_matrix PAM30 \
  --blast_word_size 2 \
  --blast_gapopen 9 \
  --blast_gapextend 1 \
  --tblastn_threshold 16 \
  --tblastn_comp_based 0 \
  --tblastn_seg no \
  -resume
```

This reproduces the equivalent of:
```
tblastn -query peptide.fasta -db nucl_db \
  -matrix PAM30 -word_size 2 \
  -gapopen 9 -gapextend 1 \
  -threshold 16 -evalue 200000 \
  -comp_based_stats 0 -seg no
```

### Step Modes

Run a specific stage of the pipeline instead of the full workflow. Useful for re-running parts with different parameters or resuming from intermediate results.

| Flag | Description |
|---|---|
| `--blast` | Run BLAST search + parse + filter only |
| `--filter` | Re-filter existing BLAST results (requires prior `--blast` run) |
| `--orf` | Run ORF finding + translation on existing filtered results. With `--mode dna`, skips codon checking and translation and writes hits directly to `standard_dna_seq/` |
| `--visual` | Run visualization (MAFFT, WebLogo, alignment figures) on existing allele results |
| `--neighbour` | Extract flanking sequence from existing `standard_dna_seq/` (requires `--upstream` and/or `--downstream`) |
| `--step <name>` | Alternative syntax: `blast`, `filter`, `orf`, `visual`, `neighbour` |

Examples:

```bash
# Run only BLAST search
nextflow run main.nf --blast --query_dir ./queries --genome_dir ./genomes --outdir ./results

# Re-filter with stricter thresholds (reuses existing BLAST .out files)
nextflow run main.nf --filter --query_dir ./queries --outdir ./results --sequence_identity 95

# Re-run visualization with different WebLogo settings
nextflow run main.nf --visual --outdir ./results --weblogo_format png --weblogo_units bits

```

**Note:** Step modes require that previous stages have already been run with the same `--outdir`. Each step reads intermediate files from the output directory.

### DNA Mode

By default the pipeline translates BLAST hits to protein before allele analysis (`--mode protein`). Set `--mode dna` to skip ORF finding and translation and work at the nucleotide level instead.

**When to use `--mode dna`:**
- Your query sequences are non-coding DNA (e.g. promoter regions, intergenic elements)
- Your sequences are pre-validated ORFs that do not need start/stop codon checking
- You want to group and compare alleles as DNA sequences rather than proteins

**What changes in `dna` mode:**
- Steps 3–4 (ORF finding + translation) are skipped entirely
- BLAST hits are used directly as allele sequences (no codon validation)
- Variant positions are called as nucleotide substitutions (e.g. `90:A→T`)
- Alignment figures use nucleotide colours (A = green, T = red, G = yellow, C = blue)
- WebLogo automatically uses the `nucleotide` color scheme

**Example:**
```bash
nextflow run main.nf \
  --query_dir ./queries \
  --genome_dir ./genomes \
  --outdir ./results \
  --mode dna
```

> **Note:** `--genome_dir` is still required in `dna` mode when `--blast_type tblastn` is used (needed to extract DNA coordinates). It is not required for `blastn` + `--mode dna`.

### WebLogo Parameters

| Parameter | Default | Description |
|---|---|---|
| `--weblogo_format` | `eps` | Output format: `eps`, `png`, `pdf`, `svg`, `jpeg` |
| `--weblogo_units` | `probability` | Y-axis units: `probability`, `bits`, `nats` |
| `--weblogo_color_scheme` | `chemistry` | Color scheme: `chemistry`, `hydrophobicity`, `charge`, `classic`, `monochrome` |
| `--weblogo_size` | `medium` | Logo size: `small`, `medium`, `large` |
| `--weblogo_resolution` | `96` | Resolution in DPI (for bitmap formats) |



## Output

Results are written to `--outdir` (default `./results/`):

```
results/
├── blast_summary/              # Raw BLAST results (CSV)
├── blast_summary_filtered/     # Filtered BLAST results (CSV)
├── standard_dna_seq/           # DNA sequences with valid start/stop codons
├── protein_seq/
│   ├── inframe/                # Longest ORF protein sequences
│   └── others/                 # Alternative reading frames (for blastn)
├── allele_info/                # KEY OUTPUT: per-query allele variant analysis
│   └── <query_name>/
│       ├── consensus_reference/  # Variants compared to the majority-vote consensus sequence
│       │   ├── allele_info_with_strain_name.csv
│       │   ├── allele_info_no_strain_name.csv
│       │   ├── allele_sequence_output.fasta
│       │   └── <query_name>_consensus.fasta
│       └── self_reference/     # Variants compared to the most abundant allele
│           ├── allele_info_with_strain_name.csv
│           ├── allele_info_no_strain_name.csv
│           └── allele_sequence_output.fasta
├── alignment_figures/           # Coloured variant alignment figures (PDF)
│   └── <query_name>/
│       ├── consensus_reference/    # Figures using consensus sequence as reference (always generated)
│       └── self_reference/         # Figures using most abundant allele as reference (only with --self_reference_figures)
├── allele_seq_filtered/        # Alleles passing species threshold
├── allele_seq/                 # MAFFT-aligned allele sequences
├── logo/                       # WebLogo sequence logos
├── blast_raw/                  # Raw BLAST output files (for --filter mode)
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
| `variant_position_count` | Number of positions differing from reference |
| `variant_percentage` | Percentage of positions that are variants |
| `variant_positions` | Specific changes, e.g. `90:K→T; 93:Q→K` (protein) or `90:A→T` (DNA) |
| `allele_sequence` | Full sequence of the allele (protein or DNA depending on `--mode`) |
| `strain_names` | Strain list (only in `_with_strain_name.csv`) |

- **`consensus_reference/`** uses a majority-vote consensus sequence (computed from all alleles) as the reference for variant calling. The CSV contains one row per allele **plus** an additional `<query_name>_consensus` reference row (with empty `species_count`). The consensus row is labeled `<query_name>_consensus` in alignment figures. Alignment figures are always generated.
- **`self_reference/`** uses the most abundant allele found as the reference. Alignment figures are generated only when `--self_reference_figures` is passed.
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

