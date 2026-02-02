
# step 4
# parse the blast results from csv to fasta
# csv_to_fasta.py first parameter input csv file, second parameter output fasta file
echo "---------------------------------------------Parsing the blast results from csv to fasta---------------------------------------------"
mkdir -p output/blast_dna_fasta
for file in output/blast_summary_filtered/*.csv; do
    echo "Processing $file"
    base=$(basename "$file" .csv)
    python code/csv_to_fasta.py "$file" "output/blast_dna_fasta/${base}.fasta"
done

################################################################################################################
# function to check start and stop codon with specified condon length
# Need two parameters
# First value as integer, for codon lenth, default to 10 codons if not specified
# Second Value as string
function check_inframe_codons() {
    local codon_range=${1:-10}
    local input_dir="$2"
    local direction=${3:-"down"}
    
    echo "---------------------------------------------Checking for inframe start codon and stop codon---------------------------------------------"
    local output_dir="output/dna_start_stop_sort_${codon_range}codon"
    mkdir -p "$output_dir"

    for file in "$input_dir"/*.fasta; do
        if [[ -f "$file" ]]; then
            echo "Processing $file"
            local base=$(basename "$file" .fasta)
            local outdir="$output_dir/${base}"
            mkdir -p "$outdir"

            python code/start_stop_codon_dna_seq.py \
                "$file" \
                "$outdir" \
                --codon_range "$codon_range" \
                --start_search_direction "$direction"
        else
            echo "No fasta files found in $input_dir/"
            break  # Exit the loop if no files found
        fi
    done
}

# function to get extended dna sequence
# need two input needed
# the first input, integer, should be the extension length in number of codon
# the second input, string, should be input dirctory (output directory of check_inframe_codons())
# the third input, string, should be output directory (replace xx with your number: output/dna_extended_xxcodon)
function get_upstream_downstream() {
    local extension_codon="$1"
    local input_dir="$2"
    local output_dir="$3"
    
    echo "---------------------------------------------Getting downstream and upstream dna sequence---------------------------------------------"
    rm -rf "$output_dir" && mkdir -p "$output_dir"

    # Loop through each subdirectory
    for dir in "$input_dir"/*/; do
        echo "Processing directory: $dir"

        dir_base=$(basename "$dir")
        output_file="$output_dir/${dir_base}.fasta"
        > "$output_file"
        
        for file in "$dir"*.fasta; do
            if [ ! -f "$file" ]; then
                continue
            fi
            
            base=$(basename "$file" .fasta)
            outfasta="$output_dir/${base}_temp.fasta"
            
            filename=$(basename "$file")
            run_mode=""

            if [[ "$filename" == "no_stop.fasta" ]]; then
                run_mode="down_stream"
            elif [[ "$filename" == "no_start.fasta" ]]; then
                run_mode="up_stream"
            elif [[ "$filename" == "no_stop_no_start.fasta" ]]; then
                run_mode="up_down_stream"
            else
                continue
            fi
            
            echo "  Processing file: $file"
            # Run the Python script
            python code/up_down_stream.py \
                --mode "$run_mode" \
                --seq_file "$file" \
                --genome_dir data/GoldenSet_7548Genomes_2023 \
                --length "$((extension_codon*3))" \
                --output_fasta "$outfasta"

            if [ -f "$outfasta" ]; then
                cat "$outfasta" >> "$output_file"
                rm "$outfasta"
            else
                echo "    Warning: Python script did not produce output for $file"
            fi
        done
    done
    
    echo "Getting upstream and downstream for no start codon and no stop codon complete"
}

# Function get the standard dna seq into a collection
# One input as strings for input file location

function get_standard_allie() {
    local input_dir="$1"
    local output_dir="output/standard_dna_seq"  # Defined output directory
    
    echo "Grabbing standard output sequences"
    mkdir -p "$output_dir"

    for dir in "$input_dir"/*/; do
        [ ! -d "$dir" ] && continue
        
        local dir_name=$(basename "$dir")        
        local source_file="$dir/standard.fasta"
        local target_file="$output_dir/${dir_name}.fasta"

        # Create empty file if target doesn't exist
        if [ ! -f "$target_file" ]; then
            touch "$target_file"
        fi

        # Check if source file exists
        if [ -f "$source_file" ]; then
            cat "$source_file" >> "$target_file" 
        else
            echo "No standard.fasta found in $dir"
        fi
    done
    
    echo "Standard sequence collection complete"
}


# step 5
# check for start and stop codon in 10 codon
    check_inframe_codons 10 output/blast_dna_fasta

# step 6
# extend the sequence by 10 codon, and check for it's stop start codon
    get_upstream_downstream 10 "output/dna_start_stop_sort_10codon" "output/dna_extended_10codon"
    check_inframe_codons 20 "output/dna_extended_10codon" "up"

# step 7
# get the standard DNA sequences
    rm -rf output/standard_dna_seq
    get_standard_allie "output/dna_start_stop_sort_10codon"
    get_standard_allie "output/dna_start_stop_sort_20codon"





# step 8
# translate standar DNA sequences into protein sequence
input_dir="output/standard_dna_seq"
out_dir="output/protein_seq"
out_inframe_dir="${out_dir}/inframe"
out_others_dir="${out_dir}/others"
mkdir -p "$out_dir"
mkdir -p "$out_inframe_dir"
mkdir -p "$out_others_dir"

echo "Processing $subdir_count subdirectories from $input_dir"

# Process each subdirectory
for file in "$input_dir"/*.fasta; do
    file_name=$(basename "$file")    
    
    # Run translation script
    python code/translate.py \
        "${input_dir}/${file_name}" \
        "${out_inframe_dir}/${file_name}" \
        "${out_others_dir}/${file_name}"
done

echo "Processing complete. $processed/$subdir_count directories processed."

# step 9
# find allels and compare it to reference query sequence
input_dir="data/query_sequences"
input_dir2="output/protein_seq/inframe"
output_dir="output/allele_info"
mkdir -p "$output_dir"

for file in "$input_dir"/*.fa*; do
    file_name=$(basename "$file" .fasta)

    # Run allel analysis script
    for infile in "${input_dir2}/${file_name}"*.fasta; do
        [[ -e "$infile" ]] || continue  # Skip if no matches
        
        # query reference compare
        out_dir_sub="${output_dir}/${file_name}/query_reference"
        mkdir -p $out_dir_sub
        python code/protein_allele_variant.py \
            -r "${input_dir}/${file_name}.fasta" \
            "$infile" \
            -o "$out_dir_sub" \
            -s "$file_name"

        # self reference compare
        out_dir_sub="${output_dir}/${file_name}/self_reference"
        mkdir -p $out_dir_sub
        python code/protein_allele_variant.py \
            "$infile" \
            -o "$out_dir_sub" \
            -s "$file_name"

    done
done

echo "Processing complete. ${output_dir} directories processed."

