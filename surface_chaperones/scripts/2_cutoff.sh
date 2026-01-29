#!/bin/bash

# Preping the enviroment 
#   To replicate the results, download code and data directory into the same directoroy
#       place the query sequence you want to check into data/quesry_sequence as .fasta file
#       add the genome sequences (e.x. files in Golden Set) into genome_sequence, also .fasta format

conda activate allele_search
cd /home/xinruoz/Document/Hiller/surface_chaperones


# step 3
# parse the blast results to csv
# filter for Evalue smaller than e-5, similarity over 80%, length 80% of query
echo "---------------------------------------------Parsing the blast results to csv---------------------------------------------"
mkdir -p output/blast_summary_filtered
input="output/blast"
input_que="data/query_sequences"
length_percentage=80
output="output/blast_summary_filtered"

for dir in "$input"/*; do
    if [ -d "$dir" ]; then
        base=$(basename "$dir")
        query="${input_que}/${base}.fasta"

        if [ ! -f "$query" ]; then
            echo "Query file not found: $query"
            continue
        fi

        # Get the length of the only sequence in the FASTA file (ignores header and line breaks)
        query_len=$(awk 'BEGIN{seq=""} !/^>/ {seq=seq$0} END{print length(seq)}' "$query")

        # Calculate 80% of the query length (rounded down)
        min_len=$((query_len * ${length_percentage} / 100))

        echo "Processing $base: query_len=$query_len, min_len=$min_len"

        python code/blastParser.py \
            --directory "$dir" \
            --querys 1 \
            --blast_type blastn \
            --max_Eval 0.00001 \
            --similarity 80 \
            --min_length "$min_len" \
            --file_name "${output}/${base}"
    fi
done


