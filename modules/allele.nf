process ALLELE_QUERY_REF {
    tag "${query_name}"
    publishDir "${params.outdir}/allele_info/${query_name}/query_reference", mode: 'copy'

    input:
    tuple val(query_name), path(protein_fasta), path(query_fasta, stageAs: 'query_ref.fasta')

    output:
    tuple val(query_name), path("*.csv"), path("allele_sequence_output.fasta")

    script:
    """
    python ${params.code_dir}/protein_allele_variant.py \
        -r query_ref.fasta \
        ${protein_fasta} \
        -o . \
        -s ${query_name}
    """
}

process ALLELE_SELF_REF {
    tag "${query_name}"
    publishDir "${params.outdir}/allele_info/${query_name}/self_reference", mode: 'copy'

    input:
    tuple val(query_name), path(protein_fasta)

    output:
    tuple val(query_name), path("*.csv"), path("allele_sequence_output.fasta"), emit: results

    script:
    """
    python ${params.code_dir}/protein_allele_variant.py \
        ${protein_fasta} \
        -o . \
        -s ${query_name}
    """
}
