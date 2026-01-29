#!/bin/bash
# Compare allele_info_with_strain_name.csv between Nextflow replication results
# and original output, listing all differences.

NF_DIR="/home/xinruoz/Document/Hiller/surface_chaperones_rep_result/results/allele_info"
OR_DIR="/home/xinruoz/Document/Hiller/surface_chaperones/output/allele_info"

for or_query_dir in "$OR_DIR"/*/; do
    query_name=$(basename "$or_query_dir")
    nf_file="$NF_DIR/$query_name/query_reference/allele_info_with_strain_name.csv"
    or_file="$OR_DIR/$query_name/query_reference/allele_info_with_strain_name.csv"

    echo "========== $query_name =========="

    if [[ ! -f "$or_file" ]]; then
        echo "  MISSING in original: $or_file"
        continue
    fi
    if [[ ! -f "$nf_file" ]]; then
        echo "  MISSING in nextflow results: $nf_file"
        continue
    fi

    # Sort both files (skip header) and diff
    header=$(head -1 "$or_file")
    diff_output=$(diff \
        <(tail -n +2 "$or_file" | sort) \
        <(tail -n +2 "$nf_file" | sort))

    if [[ -z "$diff_output" ]]; then
        echo "  IDENTICAL"
    else
        echo "  DIFFERENCES FOUND"
        # Write diff CSV: status (ONLY_IN_ORIGINAL / ONLY_IN_NEXTFLOW), then the row
        out_file="$NF_DIR/$query_name/query_reference/allele_info_with_strain_name_diff.csv"
        echo "status,$header" > "$out_file"

        # Lines only in original (< prefix)
        echo "$diff_output" | grep '^< ' | sed 's/^< //' | while IFS= read -r line; do
            echo "ONLY_IN_ORIGINAL,$line" >> "$out_file"
        done

        # Lines only in nextflow (> prefix)
        echo "$diff_output" | grep '^> ' | sed 's/^> //' | while IFS= read -r line; do
            echo "ONLY_IN_NEXTFLOW,$line" >> "$out_file"
        done

        only_orig=$(echo "$diff_output" | grep -c '^< ')
        only_nf=$(echo "$diff_output" | grep -c '^> ')
        echo "  Rows only in original: $only_orig"
        echo "  Rows only in nextflow: $only_nf"
        echo "  Diff written to: $out_file"
    fi
    echo
done
