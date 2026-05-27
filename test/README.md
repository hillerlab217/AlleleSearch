# Pipeline Test Suite

Tests for each step mode (`--filter`, `--visual`, `--blast`, `--orf`) and the full pipeline.

## Running Tests

```bash
# Tests 1-2 only (no genome data needed):
bash test/run_tests.sh

# All 5 tests (requires genome FASTA files):
bash test/run_tests.sh --genome-dir /path/to/genome/fastas
```

## What Each Test Does

| # | Mode | Needs genomes | What it checks |
|---|------|--------------|----------------|
| 1 | `--filter` | No | Re-filters pre-baked BLAST CSV; asserts filtered CSV has data |
| 2 | `--visual` | No | Generates figures/logos from pre-baked allele_info; asserts PDF + EPS exist |
| 3 | `--blast` | Yes | Full BLAST search; asserts filtered CSV has data |
| 4 | `--orf` | Yes | ORF finding on pre-baked filtered CSV; asserts DNA + protein FASTAs exist |
| 5 | full | Yes | End-to-end run; asserts all key output files exist |

Tests 1-2 use pre-baked fixture files from `test/fixtures/` (derived from the `murM/` example run).
Tests 3-5 require real genome FASTA files because BLAST and sequence extension need actual genome sequences.

## Genome Requirements (for tests 3-5)

The genome directory must contain:
- FASTA files (`.fasta`, `.fa`, `.fna`, or `.ffn`) for each genome
- File names must follow the convention `{genome_name}.fasta` — the pipeline parses species names by splitting BLAST IDs on `vs_`
- At minimum, genomes that contain the murM query locus (used in test fixtures)
- BLAST DB files (`.nin`, `.nsq`, etc.) — or pass `--make_blastdb true` to build them

The fixtures in `test/fixtures/murM_blast_filtered.csv` list the genome IDs used in the murM analysis. The genome FASTA files matching those IDs must be present in the genome directory for `--orf` (test 4) to extend sequences correctly.

## Fixtures

Pre-baked intermediate files stored in `test/fixtures/`:

| File | Source | Used by |
|------|--------|---------|
| `murM_blast_summary.csv` | `murM/blast_summary/murM.csv` | test 1 (--filter) |
| `murM_blast_filtered.csv` | `murM/blast_summary_filtered/murM.csv` | test 4 (--orf) |
| `murM_filtered.fasta` | `murM/allele_seq_filtered/murM_filtered.fasta` | test 2 (--visual) |
| `allele_info/murM/` | `murM/allele_info/murM/` | test 2 (--visual) |
