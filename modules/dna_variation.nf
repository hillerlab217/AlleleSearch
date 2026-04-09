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
