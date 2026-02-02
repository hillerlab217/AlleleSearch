process TRANSLATE {
    tag "${query_name}"
    publishDir "${params.outdir}/protein_seq", mode: 'copy'

    input:
    tuple val(query_name), path(dna_fasta)

    output:
    tuple val(query_name), path("inframe/${query_name}.fasta"), emit: inframe
    tuple val(query_name), path("others/${query_name}.fasta"), emit: others

    script:
    """
    mkdir -p inframe others
    python ${params.code_dir}/translate.py \
        ${dna_fasta} \
        inframe/${query_name}.fasta \
        others/${query_name}.fasta
    """
}
