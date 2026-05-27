process FILTER_BLAST {
    tag "${query_name}"
    publishDir "${params.outdir}/blast_summary_filtered", mode: 'copy'

    input:
    tuple val(query_name), path('input.csv'), path(query_fasta)

    output:
    tuple val(query_name), path("${query_name}.csv")

    script:
    """
    python ${params.code_dir}/filter_blast.py \
        --input input.csv \
        --output ${query_name}.csv \
        --max_evalue ${params.max_evalue} \
        --sequence_identity ${params.sequence_identity} \
        --query_fasta ${query_fasta} \
        --length_percentage ${params.length_percentage}
    """
}
