
# do weblogo and t_coffee

input="output/allele_info"
out_align="output/allele_seq"
out_logo="output/logo"
out_filtered="output/allele_seq_filtered"
out_expanded="output/allele_seq_expanded"
mkdir -p "$out_align" "$out_logo" "$out_filtered" "$out_expanded"
shopt -s nullglob

# Parameters
threshold=1        # keep sequences with species_percentage >= 1

# filter by percentage
filter_sequences_by_percentage() {
    local infile=$1  outfile=$2  thresh=$3
    awk -v threshold="$thresh" '
        /^>/ {
            keep=0
            if (match($0,/species_percentage=([0-9.]+)/,m))
                if (m[1] >= threshold) keep=1
        }
        { if (keep) print }
    ' "$infile" > "$outfile"
}

# duplicate by species count
# ---------------- Helper: duplicate by species weight (multi-line safe) -----
# • Looks for species_number=N   --> duplicates N times
# • Else falls back to rounded species_percentage
# • Keeps '*' (stop codons) or any other residue as-is
# • Rewraps sequence to 60-char lines for readability
duplicate_by_species_weight() {
    local infile=$1  outfile=$2

    awk -v wrap=60 -v outfile="$outfile" '
        BEGIN { RS=">"; ORS=""; }

        NR > 1 {
            n = split($0, L, "\n")
            header = L[1]
            seq = ""
            for (i = 2; i <= n; i++) seq = seq L[i]

            dup = 1
            if (match(header, /species_number=([0-9]+)/, m)) dup = m[1] + 0
            else if (match(header, /species_percentage=([0-9.]+)/, m))
                dup = int(m[1] + 0.5)

            for (d = 1; d <= dup; d++) {
                printf(">%s_%d\n", header, d) >> outfile
                for (p = 1; p <= length(seq); p += wrap)
                    printf("%s\n", substr(seq, p, wrap)) >> outfile
            }
        }
    ' "$infile"
}

# ---------------- Main loop -------------------
for strain_dir in "$input"/*/; do
    strain=$(basename "$strain_dir")

    if [[ -d "$strain_dir/self_reference" ]]; then
        for file in "$strain_dir/self_reference"/*.fasta; do
            echo "► $strain : $(basename "$file")"

            base="${strain}_$(basename "$file" .fasta)"
            f_filtered="${out_filtered}/${base}_filtered.fasta"
            f_aligned="${out_align}/${base}_aligned.fasta"
            f_expanded="${out_expanded}/${base}_expanded.fasta"
            logo_title="${out_logo}/${base}_logo"

            # Filter by percentage
            filter_sequences_by_percentage "$file" "$f_filtered" "$threshold"
            seq_cnt=$(grep -c '^>' "$f_filtered") || true
            if (( seq_cnt < 2 )); then
                echo "   ✘ Skipped – only $seq_cnt sequence(s) after filtering"
                continue
            fi

            # Align filtered sequences
            mafft --auto "$f_filtered" > "$f_aligned"

            # Duplicate sequences for species weighting
            duplicate_by_species_weight "$f_aligned" "$f_expanded"

            # Probability-based WebLogo
            weblogo \
                --fin  "$f_expanded" \
                --fout "$logo_title".svg \
                --format svg \
                --datatype fasta \
                --title "$base" \
                --units probability \
                --color-scheme chemistry 
        done
    fi
done



#t coffee
input_filtered="output/allele_seq_filtered"   # where filtered FASTA files live
out_align="output/allele_tcoffee"                # destination for alignments
mkdir -p "$out_align"
shopt -s nullglob                            # ignore if no matches

# ─────────────── main loop ───────────────
for f_filtered in "$input_filtered"/*_filtered.fasta; do
    base=$(basename "$f_filtered" _allele_sequence_output_filtered.fasta)
    
    out_align_sub="${out_align}/${base}"
    mkdir -p "$out_align_sub"
    # quick sanity: skip if <2 sequences
    if (( $(grep -c '^>' "$f_filtered" || true) < 2 )); then
        echo "✘ $base — skipped (fewer than 2 sequences)"
        continue
    fi

    echo "Aligning with T-Coffee: $f_filtered"
    f_clustal="${out_align_sub}/${base}_aligned_tcoffee.aln"
    f_fasta="${out_align_sub}/${base}_aligned_tcoffee.fasta"

    # T-Coffee → CLUSTAL + FASTA
    t_coffee "$f_filtered" \
        -output clustalw_aln,fasta_aln \
        -outfile "$f_clustal"

    # T-Coffee names the FASTA file *.fasta_aln; move it beside the .aln
    mv "${f_clustal%.aln}.fasta_aln" "$f_fasta"
done
