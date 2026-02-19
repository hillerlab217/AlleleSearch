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

process ALLELE_CONSENSUS_REF {
    tag "${query_name}"
    publishDir "${params.outdir}/allele_info/${query_name}/consensus_reference", mode: 'copy'

    input:
    tuple val(query_name), path(protein_fasta), path(consensus_fasta, stageAs: 'consensus_ref.fasta')

    output:
    tuple val(query_name), path("*.csv"), path("allele_sequence_output.fasta"), path("${query_name}_consensus.fasta"), emit: results

    script:
    """
    python ${params.code_dir}/protein_allele_variant.py \
        -r consensus_ref.fasta \
        ${protein_fasta} \
        -o . \
        -s ${query_name}

    cp consensus_ref.fasta ${query_name}_consensus.fasta
    """
}
