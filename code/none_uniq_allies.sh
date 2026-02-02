#!/bin/bash

# Usage: ./count_ids.sh /path/to/csv_directory output.txt

csv_dir="$1"
output_file="$2"

# Clear output file
echo "File Name, Total IDs, Unique IDs, Duplicate IDs" > "$output_file"

# Process each CSV file
for csv in "$csv_dir"/*.csv; do
    if [ ! -f "$csv" ]; then
        continue
    fi

    filename=$(basename "$csv")

    # Extract all IDs (skip header)
    ids=$(tail -n +2 "$csv" | cut -d',' -f1)

    # Count total, unique, and duplicate IDs
    total=$(echo "$ids" | wc -l)
    unique=$(echo "$ids" | sort | uniq | wc -l)
    dup=$((total - unique))

    # Print summary
    echo "$filename, $total, $unique, $dup" >> "$output_file"

    # If duplicates exist, write the duplicate lines
    if [ "$dup" -gt 0 ]; then
        echo "------ Duplicate lines in $filename ------" >> "$output_file"
        dup_ids=$(echo "$ids" | sort | uniq -d)

        # For each duplicated ID, find the lines containing it
        while IFS= read -r id; do
            grep "^$id," "$csv" >> "$output_file"
        done <<< "$dup_ids"
        echo "" >> "$output_file"
    fi
done

echo "Done. Output saved to $output_file"
