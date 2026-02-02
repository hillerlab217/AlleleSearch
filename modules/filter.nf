process FILTER_BLAST {
    tag "${query_name}"
    publishDir "${params.outdir}/blast_summary_filtered", mode: 'copy'

    input:
    tuple val(query_name), path(blast_files), path(query_fasta)

    output:
    tuple val(query_name), path("${query_name}.csv")

    script:
    """
    mkdir blast_dir
    mv ${blast_files} blast_dir/ 2>/dev/null || true

    # Calculate min_length as 80% of query sequence length
    query_len=\$(awk 'BEGIN{seq=""} !/^>/ {seq=seq\$0} END{print length(seq)}' ${query_fasta})
    min_len=\$((query_len * ${params.length_percentage} / 100))

    python ${params.code_dir}/blastParser.py \
        --directory blast_dir \
        --querys 1 \
        --blast_type ${params.blast_type ?: 'blastn'} \
        --max_Eval ${params.max_evalue} \
        --similarity ${params.similarity} \
        --min_length \$min_len \
        --file_name ${query_name}
    """
}
