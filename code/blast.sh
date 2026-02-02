#!/bin/bash

# change line 51 and line 55 if you have less than 16 cores.

# Check if correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <database_dir> <query_file_or_dir> <output_dir> <blast_type>"
    echo "  <blast_type> should be 'blastn' or 'tblastn'"
    echo "  <query_file_or_dir> can be a single .fasta file or a directory of .fasta files"
    exit 1
fi

database_dir="$1"
query_input="$2"
output_dir="$3"
blast_type="$4"

# Validate blast type
if [ "$blast_type" != "blastn" ] && [ "$blast_type" != "tblastn" ]; then
    echo "Error: Invalid blast type. Use 'blastn' or 'tblastn'."
    exit 1
fi

# Check if database directory exists
if [ ! -d "$database_dir" ]; then
    echo "Error: Database directory not found."
    exit 1
fi

# Check if query input exists (file or directory)
if [ ! -f "$query_input" ] && [ ! -d "$query_input" ]; then
    echo "Error: Query file or directory not found."
    exit 1
fi

# Create main output directory if it doesn't exist
mkdir -p "$output_dir"

# Function to run BLAST for a single query against all DB sequences
run_blast() {
    local query_file="$1"
    local query_name=$(basename "$query_file" .fasta)
    local query_output_dir="$output_dir/$query_name"
    
    mkdir -p "$query_output_dir"

    # Process each database file in the database directory
    for db_file in "$database_dir"/*.{fasta,fa,fna,ffn}; do
        [ -f "$db_file" ] || continue  # Skip if no matching files
        
        local db_name=$(basename "$db_file")
        local db_prefix="${db_name%.*}"
        local out_file="$query_output_dir/${query_name}_vs_${db_prefix}.out"
        
        shopt -s nullglob
        db_prefix="${db_file%.*}"

        # Rebuild BLAST database if any of the required index files are missing or empty
        if [[ ! -s "${db_prefix}.nhr" || ! -s "${db_prefix}.nin" || ! -s "${db_prefix}.nsq" ]]; then
            echo "Building BLAST DB for $db_file ..."
            rm -f "${db_prefix}".{nhr,nin,nsq} 2>/dev/null
            makeblastdb -in "$db_file" -dbtype nucl -out "$db_prefix" -parse_seqids
        fi

        echo "Processing $query_name against $(basename "$db_prefix")..."

        if [ "$blast_type" == "blastn" ]; then
            blastn -query "$query_file" \
                -db "$db_prefix" \
                -out "$out_file" \
                -evalue 1e-5 \
                -num_threads 4 \
                -outfmt 6
        else
            tblastn -query "$query_file" \
                    -db "$db_prefix" \
                    -out "$out_file" \
                    -evalue 1e-5 \
                    -num_threads 4 \
                    -outfmt 6
        fi

    done
}

# If input is a single file, process it
if [ -f "$query_input" ]; then
    run_blast "$query_input"
# If input is a directory, process all .fasta files inside it
elif [ -d "$query_input" ]; then
    for query_file in "$query_input"/*.{fasta,fa,fna,ffn}; do
        [ -f "$query_file" ] || continue  # Skip if no matching files
        echo "Processing $query_file..."
        run_blast "$query_file"
    done
fi

echo "${blast_type^} execution completed. Results saved in $output_dir/"