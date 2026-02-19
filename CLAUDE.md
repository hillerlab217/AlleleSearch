# Allele Search Pipeline — CLAUDE.md

## Project Purpose

This is a **Nextflow DSL2 bioinformatics pipeline** for discovering protein allele variants across large collections of bacterial/microbial genomes. Given one or more query DNA/protein sequences and a directory of genome FASTA files, it:

1. BLASTs each query against all genomes
2. Validates and extends open reading frames (start/stop codon checking)
3. Translates DNA hits to proteins
4. Clusters identical proteins into allele groups
5. Calls amino acid variants relative to a reference
6. Produces colored alignment figures and WebLogo sequence logos

**Domain context:** "Allele" here means a unique protein variant found in ≥1 sequenced strain. Variant positions are formatted as `90:K→T` (position : reference AA → observed AA). The most common allele is `allele_1`.

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

The `nextflow.config` hard-codes the conda env path:
```
export PATH=/home/xinruoz/miniconda3/envs/allele_search/bin:$PATH
```
IMPORTANT: Update this path in `nextflow.config` if the conda prefix changes.

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

### Step-specific modes (partial runs)
Each mode reads from and writes to `--outdir`. Prior steps must have been run first.
```bash
# Re-filter with different thresholds (reuses existing BLAST .out files)
nextflow run main.nf --filter --query_dir ./data --outdir ./results --sequence_identity 95

# Re-run ORF finding + translation
nextflow run main.nf --orf --query_dir ./data --genome_dir /genomes --outdir ./results

# Re-run visualization only
nextflow run main.nf --visual --outdir ./results --weblogo_format png
```

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
main.nf                     # Workflow entrypoint; routes --step modes
nextflow.config             # All parameters with defaults
environment.yml             # Conda deps
modules/
  blast.nf                  # COLLECT_ALL_BLAST_DBS, BLAST_SEARCH, PARSE_BLAST
  filter.nf                 # FILTER_BLAST
  orf.nf                    # EXTRACT_DNA, CSV_TO_FASTA, CODON_CHECK, EXTEND_SEQUENCES, COLLECT_STANDARD
  translate.nf              # TRANSLATE (DNA → protein, 3 frames)
  allele.nf                 # ALLELE_QUERY_REF, ALLELE_SELF_REF
  quality.nf                # FILTER_ALLELES, MAFFT_ALIGN, WEBLOGO, ALIGNMENT_FIGURE
code/
  blastParser.py            # Converts raw BLAST .out → CSV
  csv_to_fasta.py           # CSV → FASTA with metadata headers
  start_stop_codon_dna_seq.py  # Classifies sequences by codon validity
  up_down_stream.py         # Extends incomplete sequences from genome files
  translate.py              # DNA → protein in 3 frames, returns longest ORF
  extract_dna_from_blast_csv.py  # Replaces tblastn protein seqs with DNA coords
  protein_allele_variant.py # Core: groups identical proteins, calls variants (BLOSUM62)
  alignment_figure.py       # Renders colored AA alignment PDFs/PNGs
data/                       # Example query FASTA files (insertion regions)
test_query/                 # Test data (murM.fa)
results/                    # Pipeline outputs (gitignored)
```

---

## Key Conventions & Design Decisions

### FASTA header format (internal)
Scripts communicate metadata through structured FASTA headers. DO NOT change this format without updating all downstream parsers:
```
>ID|contig=NAME|E-value=VALUE|AlignLen=LEN|SeqIdentity=%|SubjectStart=POS|SubjectEnd=POS
```

### Strand detection
**SubjectStart > SubjectEnd signals reverse strand.** This is used in `up_down_stream.py` and `extract_dna_from_blast_csv.py`. No explicit strand column exists.

### Species name extraction
Species/strain names are parsed from BLAST IDs by splitting on `"vs_"`. BLAST output file names follow the convention `{query}_vs_{genome_name}`. Changing this naming convention will break allele species counting in `protein_allele_variant.py`.

### Two allele reference modes (always both generated)
- `query_reference/`: variants called against the original query sequence
- `self_reference/`: variants called against the most abundant allele found

Both are always produced in the full pipeline. The self-reference is used for WebLogo weighting.

### BLAST batching
Genomes are batched 100 per BLAST task (configurable: `--blast_batch_size`). Each batch produces a separate `.out` file; `PARSE_BLAST` collects all batches per query with `groupTuple()`.

### Codon check logic
- Round 1 (`CODON_CHECK_10`): 10-codon window at each terminus
- Round 2 (`CODON_CHECK_20`): 20-codon window after genome extension
- Valid start codons: ATG, GTG, TTG
- Valid stop codons: TAA, TAG, TGA
- Sequences with internal in-frame stop codons are discarded

### Amino acid color scheme (alignment figures)
Physicochemical groups in `alignment_figure.py`:
- Orange `#f4b183`: Hydrophobic (A G I L M V P)
- Blue `#4f81bd`: Positive (R K H)
- Green `#9bbb59`: Polar (N C Q S T)
- Red `#e26b6b`: Negative (D E)
- Yellow `#fff475`: Aromatic (F W Y)

---

## Common Gotchas

- **`--filter` mode requires `--query_dir`** — needed to calculate the minimum alignment length from query sequence length. Forgetting this causes an error.
- **`--blast_type tblastn` triggers an extra step** (`EXTRACT_DNA`) to swap protein sequences back for DNA coordinates before ORF checking.
- **WebLogo weighting** in `quality.nf` expands sequences by `species_percentage` (rounded to int) to weight the logo by strain frequency. Headers must contain `species_percentage=X.XX`.
- **MAFFT is skipped silently** if fewer than 2 sequences pass the `--species_threshold` filter. Check `WARNING` messages in nextflow logs.
- **Nextflow work directory** can grow very large. Run `nextflow clean -f` periodically to remove cached task directories.
- **`-resume` reuses cached tasks** — omit it when you want to force a full re-run. Step modes (`--blast`, `--filter`, etc.) always overwrite their output directories regardless of `-resume`.
- **BlastDB building** is skipped by default (`--make_blastdb false`). If `.n*` index files are missing, set `--make_blastdb true`.

---

## Development Workflow

### Adding a new parameter
1. Add to `params {}` block in `nextflow.config` with a sensible default
2. Reference via `${params.param_name}` in the relevant module's `script:` block
3. Document in `README.md` parameter tables
4. Add to `--help` message in `main.nf`

### Adding a new pipeline step
1. Create or edit the relevant `modules/*.nf` file
2. Import the process in `main.nf` at the top
3. Wire it into the appropriate workflow block(s) (`full`, and any relevant step modes)
4. If it produces output needed by downstream steps, update channel routing

### Modifying Python scripts
Scripts in `code/` are called from Nextflow `script:` blocks via:
```
python ${params.code_dir}/script_name.py [args]
```
The `code_dir` defaults to `${projectDir}/code`. Scripts must be executable and accept arguments from the command line (argparse or sys.argv).

### Testing a step change
```bash
# Re-run only the affected step with -resume to reuse upstream work
nextflow run main.nf --visual --outdir ./results/test -resume
```

---

## Output Interpretation

Key output: `results/allele_info/<query_name>/`

| File | Description |
|---|---|
| `self_reference/allele_info_no_strain_name.csv` | Alleles ranked by frequency; variants vs. most common allele |
| `query_reference/allele_info_no_strain_name.csv` | Same, but variants vs. original query sequence |
| `allele_sequence_output.fasta` | Protein FASTA, one entry per allele, header includes species count |

CSV columns: `allele_name`, `species_count`, `species_percentage`, `variant_position_count`, `variant_percentage`, `variant_positions(position:reference→changed)`, `allele_sequence`

Alignment figures: `results/alignment_figures/<query>/{self,query}_reference/`
- `*_alignment.pdf` — only variant positions shown
- `*_alignment_gap.pdf` — all positions between first and last variant

Logos: `results/logo/<query_name>_logo.<format>`

---

## Related Documentation

- `README.md` — User-facing installation and usage guide
- `run_example.txt` — Real command examples from previous analysis runs
- `task.md` — Active development task list
- `flowchart.mermaid` — Visual overview of the pipeline DAG
