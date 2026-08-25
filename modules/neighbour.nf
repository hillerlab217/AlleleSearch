process EXTRACT_NEIGHBOURS {
    tag "${query_name}"
    publishDir "${params.outdir}/neighbour/${query_name}", mode: 'copy'

    input:
    tuple val(query_name), path(dna_fasta)
    path genome_dir

    output:
    tuple val(query_name), path("*.fasta")

    script:
    def merge_flag = params.neighbour_merge ? '--merge'      : ''
    def flank_flag = params.flank_only      ? '--flank_only' : ''
    """
    python ${params.code_dir}/extract_neighbours.py \
        --seq_file ${dna_fasta} \
        --genome_dir ${genome_dir} \
        --upstream ${params.upstream} \
        --downstream ${params.downstream} \
        ${merge_flag} ${flank_flag} \
        --output_dir .
    """
}
