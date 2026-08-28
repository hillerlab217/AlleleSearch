// ── BLAST processes that always use blastn (DNA conservation mode) ─────────

process BLAST_SEARCH_DNA {
    tag "${query.baseName}_batch"
    cpus 2
    publishDir "${params.outdir}/dna_variation/blast_raw", mode: 'copy'

    input:
    tuple path(query), val(genome_names_csv), path(genome_dir)

    output:
    tuple val("${query.baseName}"), path("*.out")

    script:
    def shared_args = [
        params.blast_word_size != null ? "-word_size ${params.blast_word_size}" : '',
        params.blast_gapopen   != null ? "-gapopen ${params.blast_gapopen}"     : '',
        params.blast_gapextend != null ? "-gapextend ${params.blast_gapextend}" : '',
    ].findAll { it }.join(' ')
    """
    for gname in \$(echo '${genome_names_csv}' | tr ',' ' '); do
        genome="${genome_dir}/\${gname}"
        [ ! -e "\$genome" ] && continue
        base=\$(echo "\$gname" | sed 's/\\.[^.]*\$//')

        if [ -f "\${genome}.nin" ] || [ -f "\${genome}.nsq" ]; then
            DB_PREFIX="\${genome}"
        else
            DB_PREFIX="\${genome%.*}"
        fi

        blastn -query ${query} \
            -db \${DB_PREFIX} \
            -out ${query.baseName}_vs_\${base}.out \
            -evalue ${params.max_evalue} \
            ${shared_args} \
            -num_threads ${task.cpus}
    done
    """
}

process PARSE_BLAST_DNA {
    tag "${query_name}"
    publishDir "${params.outdir}/dna_variation/blast_summary", mode: 'copy'

    input:
    tuple val(query_name), path(blast_files)

    output:
    tuple val(query_name), path("${query_name}.csv")

    script:
    """
    mkdir blast_dir
    cp ${blast_files} blast_dir/ 2>/dev/null || true
    python ${params.code_dir}/blastParser.py \
        --directory blast_dir \
        --querys 1 \
        --blast_type blastn \
        --max_Eval 10 \
        --sequence_identity 0 \
        --min_length 0 \
        --file_name ${query_name}
    """
}

process FILTER_BLAST_DNA {
    tag "${query_name}"
    publishDir "${params.outdir}/dna_variation/blast_summary_filtered", mode: 'copy'

    input:
    tuple val(query_name), path(input_csv), path(query_fasta)

    output:
    tuple val(query_name), path("filtered.csv")

    script:
    """
    python ${params.code_dir}/filter_blast.py \
        --input ${input_csv} \
        --output filtered.csv \
        --max_evalue ${params.max_evalue} \
        --sequence_identity ${params.sequence_identity} \
        --length_percentage ${params.length_percentage} \
        --query_fasta ${query_fasta}
    """
}

// ── MAFFT + conservation figure ────────────────────────────────────────────

process MAFFT_DNA {
    tag "${query_name}"
    publishDir "${params.outdir}/dna_variation/mafft", mode: 'copy'

    input:
    tuple val(query_name), path(dna_fasta)

    output:
    tuple val(query_name), path("${query_name}_dna_mafft.fasta")

    script:
    """
    seq_cnt=\$(grep -c '^>' ${dna_fasta} || true)
    if [ "\$seq_cnt" -ge 2 ]; then
        mafft --auto ${dna_fasta} > ${query_name}_dna_mafft.fasta
    else
        echo "Skipping MAFFT: only \$seq_cnt sequence(s) in ${dna_fasta}" >&2
        cp ${dna_fasta} ${query_name}_dna_mafft.fasta
    fi
    """
}

process DNA_CONSERVATION_FIGURE {
    tag "${query_name}"
    publishDir "${params.outdir}/dna_variation/figures", mode: 'copy'

    input:
    tuple val(query_name), path(aligned_fasta)

    output:
    tuple val(query_name), path("${query_name}_conservation.${params.dna_variation_figure_format}"), optional: true

    script:
    """
    python ${params.code_dir}/dna_conservation_figure.py \
        ${aligned_fasta} \
        ${query_name}_conservation.${params.dna_variation_figure_format} \
        ${query_name} \
        --format ${params.dna_variation_figure_format} \
    || echo "Skipping conservation figure for ${query_name}" >&2
    """
}
