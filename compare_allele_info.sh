#!/bin/bash
# Compare CSV files between rep_result and output allele_info directories
# For each CSV that exists in rep_result, show differences with the corresponding CSV in output

REP_DIR="/home/xinruoz/Document/Hiller/allele_search/surface_chaperones_rep_result/allele_info"
OUT_DIR="/home/xinruoz/Document/Hiller/allele_search/surface_chaperones/output/allele_info"
RESULTS="/home/xinruoz/Document/Hiller/allele_search/compare_allele_info_results.txt"

{
for gene_dir in "$REP_DIR"/*/; do
    gene=$(basename "$gene_dir")

    # Check if corresponding gene dir exists in output
    if [ ! -d "$OUT_DIR/$gene" ]; then
        echo "=== $gene: directory missing in output ==="
        echo ""
        continue
    fi

    for ref_type in query_reference self_reference; do
        for csv in "$REP_DIR/$gene/$ref_type"/*.csv; do
            [ -f "$csv" ] || continue
            csv_name=$(basename "$csv")
            other="$OUT_DIR/$gene/$ref_type/$csv_name"

            echo "=== $gene / $ref_type / $csv_name ==="

            if [ ! -f "$other" ]; then
                echo "  MISSING in output"
            else
                diff_output=$(diff "$csv" "$other")
                if [ -z "$diff_output" ]; then
                    echo "  IDENTICAL"
                else
                    echo "$diff_output"
                fi
            fi
            echo ""
        done
    done
done
} | tee "$RESULTS"

echo ""
echo "Results written to $RESULTS"
