#!/usr/bin/env bash
# run_tests.sh — Test suite for allele_search pipeline step modes
#
# Usage:
#   bash test/run_tests.sh                          # runs filter + visual tests only
#   bash test/run_tests.sh --genome-dir /path/to/genomes   # runs all 5 tests
#
# Tests 1-2 (--filter, --visual) run without genome data using pre-baked fixtures.
# Tests 3-5 (--blast, --orf, full) require actual genome FASTA files.
# See test/README.md for genome requirements.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
FIXTURES_DIR="$SCRIPT_DIR/fixtures"
TEST_OUTDIR="$SCRIPT_DIR/test_output"
QUERY_DIR="$PIPELINE_DIR/test_query"
GENOME_DIR=""

PASS=0
FAIL=0
SKIP=0

# --- Argument parsing ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome-dir) GENOME_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# --- Helpers ---
print_result() {
    local status="$1" name="$2" msg="${3:-}"
    if   [[ "$status" == "PASS" ]]; then echo "  [PASS] $name"; ((PASS++))
    elif [[ "$status" == "FAIL" ]]; then echo "  [FAIL] $name${msg:+: $msg}"; ((FAIL++))
    elif [[ "$status" == "SKIP" ]]; then echo "  [SKIP] $name${msg:+: $msg}"; ((SKIP++))
    fi
}

assert_file_exists() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        echo "    Missing: $path"
        return 1
    fi
}

assert_csv_has_data() {
    local path="$1"
    local lines
    lines=$(wc -l < "$path")
    if [[ "$lines" -lt 2 ]]; then
        echo "    CSV has no data rows: $path"
        return 1
    fi
}

seed_dir() {
    # Copy a fixture file/dir into the test outdir, creating parent dirs as needed
    local src="$1" dest="$2"
    mkdir -p "$(dirname "$dest")"
    cp -r "$src" "$dest"
}

cleanup() {
    rm -rf "$TEST_OUTDIR"
}

# --- Tests ---

test_filter() {
    local outdir="$TEST_OUTDIR/filter"
    rm -rf "$outdir"
    seed_dir "$FIXTURES_DIR/murM_blast_summary.csv" "$outdir/blast_summary/murM.csv"

    if nextflow run "$PIPELINE_DIR/main.nf" \
        --filter \
        --query_dir "$QUERY_DIR" \
        --outdir "$outdir" \
        --sequence_identity 80 \
        --length_percentage 80 \
        -ansi-log false \
        > "$outdir/nextflow.log" 2>&1; then
        if assert_file_exists "$outdir/blast_summary_filtered/murM.csv" && \
           assert_csv_has_data  "$outdir/blast_summary_filtered/murM.csv"; then
            print_result PASS "--filter mode"
        else
            print_result FAIL "--filter mode" "output file missing or empty"
        fi
    else
        print_result FAIL "--filter mode" "nextflow exited with error (see $outdir/nextflow.log)"
    fi
}

test_visual() {
    local outdir="$TEST_OUTDIR/visual"
    rm -rf "$outdir"
    seed_dir "$FIXTURES_DIR/murM_filtered.fasta"             "$outdir/allele_seq_filtered/murM_filtered.fasta"
    seed_dir "$FIXTURES_DIR/allele_info/murM/self_reference"      "$outdir/allele_info/murM/self_reference"
    seed_dir "$FIXTURES_DIR/allele_info/murM/consensus_reference" "$outdir/allele_info/murM/consensus_reference"

    if nextflow run "$PIPELINE_DIR/main.nf" \
        --visual \
        --outdir "$outdir" \
        -ansi-log false \
        > "$outdir/nextflow.log" 2>&1; then
        local ok=1
        assert_file_exists "$outdir/logo/murM_logo.eps"                                          || ok=0
        assert_file_exists "$outdir/alignment_figures/murM/self_reference/murM_alignment.pdf"    || ok=0
        assert_file_exists "$outdir/alignment_figures/murM/consensus_reference/murM_alignment.pdf" || ok=0
        if [[ $ok -eq 1 ]]; then
            print_result PASS "--visual mode"
        else
            print_result FAIL "--visual mode" "one or more output files missing"
        fi
    else
        print_result FAIL "--visual mode" "nextflow exited with error (see $outdir/nextflow.log)"
    fi
}

test_blast() {
    local outdir="$TEST_OUTDIR/blast"
    rm -rf "$outdir"

    if nextflow run "$PIPELINE_DIR/main.nf" \
        --blast \
        --query_dir "$QUERY_DIR" \
        --genome_dir "$GENOME_DIR" \
        --outdir "$outdir" \
        --make_blastdb true \
        -ansi-log false \
        > "$outdir/nextflow.log" 2>&1; then
        if assert_file_exists "$outdir/blast_summary_filtered/murM.csv" && \
           assert_csv_has_data  "$outdir/blast_summary_filtered/murM.csv"; then
            print_result PASS "--blast mode"
        else
            print_result FAIL "--blast mode" "output file missing or empty"
        fi
    else
        print_result FAIL "--blast mode" "nextflow exited with error (see $outdir/nextflow.log)"
    fi
}

test_orf() {
    local outdir="$TEST_OUTDIR/orf"
    rm -rf "$outdir"
    seed_dir "$FIXTURES_DIR/murM_blast_filtered.csv" "$outdir/blast_summary_filtered/murM.csv"

    if nextflow run "$PIPELINE_DIR/main.nf" \
        --orf \
        --query_dir "$QUERY_DIR" \
        --genome_dir "$GENOME_DIR" \
        --outdir "$outdir" \
        -ansi-log false \
        > "$outdir/nextflow.log" 2>&1; then
        local ok=1
        assert_file_exists "$outdir/standard_dna_seq/murM.fasta"    || ok=0
        assert_file_exists "$outdir/protein_seq/inframe/murM.fasta"  || ok=0
        if [[ $ok -eq 1 ]]; then
            print_result PASS "--orf mode"
        else
            print_result FAIL "--orf mode" "one or more output files missing"
        fi
    else
        print_result FAIL "--orf mode" "nextflow exited with error (see $outdir/nextflow.log)"
    fi
}

test_full() {
    local outdir="$TEST_OUTDIR/full"
    rm -rf "$outdir"

    if nextflow run "$PIPELINE_DIR/main.nf" \
        --query_dir "$QUERY_DIR" \
        --genome_dir "$GENOME_DIR" \
        --outdir "$outdir" \
        --make_blastdb true \
        -ansi-log false \
        > "$outdir/nextflow.log" 2>&1; then
        local ok=1
        assert_file_exists "$outdir/blast_summary_filtered/murM.csv"                                          || ok=0
        assert_file_exists "$outdir/standard_dna_seq/murM.fasta"                                             || ok=0
        assert_file_exists "$outdir/protein_seq/inframe/murM.fasta"                                          || ok=0
        assert_file_exists "$outdir/allele_info/murM/self_reference/allele_info_no_strain_name.csv"          || ok=0
        assert_file_exists "$outdir/allele_info/murM/consensus_reference/allele_info_no_strain_name.csv"     || ok=0
        assert_file_exists "$outdir/allele_seq_filtered/murM_filtered.fasta"                                 || ok=0
        assert_file_exists "$outdir/logo/murM_logo.eps"                                                      || ok=0
        if [[ $ok -eq 1 ]]; then
            print_result PASS "full pipeline"
        else
            print_result FAIL "full pipeline" "one or more output files missing"
        fi
    else
        print_result FAIL "full pipeline" "nextflow exited with error (see $outdir/nextflow.log)"
    fi
}

# --- Main ---
echo "Allele Search Pipeline — Test Suite"
echo "===================================="
echo ""

echo "Tests that run without genome data:"
test_filter
test_visual

echo ""
if [[ -z "$GENOME_DIR" ]]; then
    echo "Tests requiring genome data (skipped — pass --genome-dir to enable):"
    print_result SKIP "--blast mode"  "requires --genome-dir"
    print_result SKIP "--orf mode"    "requires --genome-dir"
    print_result SKIP "full pipeline" "requires --genome-dir"
else
    echo "Tests requiring genome data (genome-dir: $GENOME_DIR):"
    test_blast
    test_orf
    test_full
fi

echo ""
echo "===================================="
echo "Results: $PASS passed, $FAIL failed, $SKIP skipped"
echo ""

if [[ $FAIL -gt 0 ]]; then
    echo "Some tests failed. Check log files in $TEST_OUTDIR/*/nextflow.log"
    exit 1
fi
