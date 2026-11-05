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
    if   [[ "$status" == "PASS" ]]; then echo "  [PASS] $name"; (( ++PASS ))
    elif [[ "$status" == "FAIL" ]]; then echo "  [FAIL] $name${msg:+: $msg}"; (( ++FAIL ))
    elif [[ "$status" == "SKIP" ]]; then echo "  [SKIP] $name${msg:+: $msg}"; (( ++SKIP ))
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
        --length_percentage 30 \
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
        --self_reference_figures true \
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

test_smoke() {
    # Full pipeline on two tiny synthetic genomes (wildtype + K86N).
    # Uses test/fixtures/genomes/ — no external genome data required.
    # Asserts:
    #   1. allele CSV produced with ≥2 alleles
    #   2. One allele has zero variant positions (wildtype)
    #   3. One allele's variant field contains "86:K→N" (the K86N genome)
    local outdir="$TEST_OUTDIR/smoke"
    local toy_genomes="$SCRIPT_DIR/fixtures/genomes"
    rm -rf "$outdir"
    mkdir -p "$outdir"

    if [[ ! -d "$toy_genomes" ]]; then
        print_result FAIL "smoke test" "fixtures/genomes/ missing — run: python test/make_fixtures.py"
        return
    fi
    if [[ ! -f "$QUERY_DIR/murM.fa" ]]; then
        print_result FAIL "smoke test" "test_query/murM.fa missing — run: python test/make_fixtures.py"
        return
    fi

    if nextflow run "$PIPELINE_DIR/main.nf" \
        --query_dir "$QUERY_DIR" \
        --genome_dir "$toy_genomes" \
        --outdir "$outdir" \
        --make_blastdb true \
        -ansi-log false \
        > "$outdir/nextflow.log" 2>&1; then

        local csv="$outdir/allele_info/murM/self_reference/allele_info_no_strain_name.csv"

        if ! assert_file_exists "$csv"; then
            print_result FAIL "smoke test" "allele CSV missing"
            return
        fi

        # Must have ≥2 allele rows
        local nrows
        nrows=$(tail -n +2 "$csv" | wc -l)
        if [[ "$nrows" -lt 2 ]]; then
            print_result FAIL "smoke test" "expected ≥2 alleles, got $nrows"
            return
        fi

        # One allele must have variant_position_count=0 (wildtype)
        local wt_rows
        wt_rows=$(awk -F',' 'NR>1 && $4=="0"' "$csv" | wc -l)
        if [[ "$wt_rows" -lt 1 ]]; then
            print_result FAIL "smoke test" "no wildtype allele (variant_position_count=0) found in CSV"
            return
        fi

        # One allele must contain "86:K→N" in its variant field (K86N genome)
        # variant_positions column is column 6 (0-indexed 5); use Python to handle quoting safely
        if python3 - "$csv" <<'PYEOF'
import csv, sys
with open(sys.argv[1]) as f:
    for row in csv.DictReader(f):
        vfield = row.get('variant_positions(position:reference\u2192changed)', '')
        if '86:K\u2192N' in vfield:
            sys.exit(0)
sys.exit(1)
PYEOF
        then
            print_result PASS "smoke test"
        else
            print_result FAIL "smoke test" "variant '86:K→N' not found in any allele row"
        fi

    else
        print_result FAIL "smoke test" "nextflow exited with error (see $outdir/nextflow.log)"
    fi
}

test_known_answer() {
    # Known-answer E2E test: planted K86N mutation must round-trip through the
    # full pipeline and appear as "86:K→N" in the consensus-reference allele CSV.
    #
    # Setup: 3 genomes:
    #   wildtype.fasta   → allele with 0 variants (reference)
    #   wildtype_v2.fasta→ same wildtype (simulates a second strain)
    #   K86N.fasta        → planted mutation at position 86 (K→N)
    #
    # With 2 wildtype and 1 K86N genome:
    #   - Consensus sequence at position 86 = K (majority)
    #   - allele_1 = wildtype (2 strains, 66.67%), 0 variants vs consensus
    #   - allele_2 = K86N    (1 strain,  33.33%), variant "86:K→N" vs consensus
    local outdir="$TEST_OUTDIR/known_answer"
    local toy_genomes="$SCRIPT_DIR/fixtures/genomes"
    local ka_genomes="$TEST_OUTDIR/known_answer_genomes"
    rm -rf "$outdir" "$ka_genomes"
    mkdir -p "$outdir" "$ka_genomes"

    if [[ ! -f "$toy_genomes/wildtype.fasta" || ! -f "$toy_genomes/K86N.fasta" ]]; then
        print_result FAIL "known-answer test" "fixtures missing — run: python test/make_fixtures.py"
        return
    fi

    # Provide exactly 2 wildtype + 1 K86N genomes for deterministic allele counts
    cp "$toy_genomes/wildtype.fasta" "$ka_genomes/wildtype.fasta"
    cp "$toy_genomes/wildtype.fasta" "$ka_genomes/wildtype_v2.fasta"
    cp "$toy_genomes/K86N.fasta"     "$ka_genomes/K86N.fasta"

    if nextflow run "$PIPELINE_DIR/main.nf" \
        --query_dir "$QUERY_DIR" \
        --genome_dir "$ka_genomes" \
        --outdir "$outdir" \
        --make_blastdb true \
        -ansi-log false \
        > "$outdir/nextflow.log" 2>&1; then

        local cons_csv="$outdir/allele_info/murM/consensus_reference/allele_info_no_strain_name.csv"
        local self_csv="$outdir/allele_info/murM/self_reference/allele_info_no_strain_name.csv"
        local ok=1

        # Both CSVs must exist
        assert_file_exists "$cons_csv" || ok=0
        assert_file_exists "$self_csv" || ok=0
        if [[ $ok -eq 0 ]]; then
            print_result FAIL "known-answer test" "allele CSV(s) missing"
            return
        fi

        # Exactly 2 allele rows (3 genomes → 2 unique alleles).
        # The consensus-reference CSV also contains a murM_consensus reference row
        # (species_count column is empty); exclude it when counting allele rows.
        local nrows
        nrows=$(awk -F',' 'NR>1 && $2!=""' "$cons_csv" | grep -v '^$' | wc -l)
        if [[ "$nrows" -ne 2 ]]; then
            print_result FAIL "known-answer test" "expected 2 allele rows in consensus CSV, got $nrows"
            return
        fi

        # allele_1 must have 0 variant positions
        local wt_rows
        wt_rows=$(awk -F',' 'NR>1 && $4=="0"' "$cons_csv" | wc -l)
        if [[ "$wt_rows" -lt 1 ]]; then
            print_result FAIL "known-answer test" "allele_1 should have 0 variants vs consensus"
            return
        fi

        # allele_2 must have "86:K→N" in consensus-reference CSV
        if ! python3 - "$cons_csv" <<'PYEOF'
import csv, sys
with open(sys.argv[1]) as f:
    for row in csv.DictReader(f):
        vfield = row.get('variant_positions(position:reference\u2192changed)', '')
        if '86:K\u2192N' in vfield:
            sys.exit(0)
sys.exit(1)
PYEOF
        then
            print_result FAIL "known-answer test" "86:K\u2192N not found in consensus-reference CSV"
            return
        fi

        # allele_2 must also have "86:K→N" in self-reference CSV
        if ! python3 - "$self_csv" <<'PYEOF'
import csv, sys
with open(sys.argv[1]) as f:
    for row in csv.DictReader(f):
        vfield = row.get('variant_positions(position:reference\u2192changed)', '')
        if '86:K\u2192N' in vfield:
            sys.exit(0)
sys.exit(1)
PYEOF
        then
            print_result FAIL "known-answer test" "86:K\u2192N not found in self-reference CSV"
            return
        fi

        print_result PASS "known-answer test"
    else
        print_result FAIL "known-answer test" "nextflow exited with error (see $outdir/nextflow.log)"
    fi
}

# --- Main ---
echo "Allele Search Pipeline — Test Suite"
echo "===================================="
echo ""

echo "Tests using bundled toy fixtures (no external genome data needed):"
test_filter
test_visual
test_smoke
test_known_answer

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
