#!/bin/bash

# Preping the enviroment 
#   To replicate the results, download code and data directory into the same directoroy
#       place the query sequence you want to check into data/quesry_sequence as .fasta file
#       add the genome sequences (e.x. files in Golden Set) into genome_sequence, also .fasta format

conda activate allele_search
cd /home/xinruoz/Document/Hiller/surface_chaperones

# step 1
# blast for the targeted proteins with DNA sequence from D39v
echo "---------------------------------------------Blasting---------------------------------------------"
mkdir -p output
mkdir -p output/blast
bash code/blast.sh data/genome_sequence query_sequences output/blast blastn

# step 2
# parse the blast results to csv without filtering
echo "---------------------------------------------Parsing the blast results to csv without filtering---------------------------------------------"
mkdir -p output/blast_summary
input="output/blast"
output="output/blast_summary"

for dir in "$input"/*; do
    base=$(basename "$dir")
  python code/blastParser.py \
    --directory "$dir" \
    --querys 1 \
    --blast_type blastn \
    --max_Eval 10 \
    --similarity 0 \
    --min_length 0 \
    --file_name "${output}/${base}"
done

bash code/none_uniq_allies.sh output/blast_summary output/blast_summary/none_uniq_allies_in_species.txt

echo "---------------------------------------------Done!---------------------------------------------"
echo "Check output/blast_summary folder for the results"