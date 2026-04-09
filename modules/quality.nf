process MAFFT_FOR_CONSENSUS {
    tag "${query_name}"
    publishDir "${params.outdir}/allele_mafft_consensus", mode: 'copy'

    input:
    tuple val(query_name), path(allele_fasta)

    output:
    tuple val(query_name), path("${query_name}_all_aligned.fasta")

    script:
    """
    seq_cnt=\$(grep -c '^>' ${allele_fasta} || true)
    if [ "\$seq_cnt" -ge 2 ]; then
        mafft --auto ${allele_fasta} > ${query_name}_all_aligned.fasta
    else
        cp ${allele_fasta} ${query_name}_all_aligned.fasta
    fi
    """
}

process COMPUTE_CONSENSUS {
    tag "${query_name}"

    input:
    tuple val(query_name), path(aligned_fasta)

    output:
    tuple val(query_name), path("${query_name}_consensus.fasta")

    script:
    """
    python ${params.code_dir}/compute_consensus.py \
        ${aligned_fasta} \
        -n ${query_name} \
        -o ${query_name}_consensus.fasta
    """
}

process FILTER_ALLELES {
    tag "${query_name}"
    publishDir "${params.outdir}/allele_seq_filtered", mode: 'copy'

    input:
    tuple val(query_name), path(allele_fasta)

    output:
    tuple val(query_name), path("${query_name}_filtered.fasta")

    script:
    """
    awk -v threshold=${params.species_threshold} '
        /^>/ {
            keep=0
            if (match(\$0,/species_percentage=([0-9.]+)/,m))
                if (m[1] >= threshold) keep=1
        }
        { if (keep) print }
    ' ${allele_fasta} > ${query_name}_filtered.fasta

    # Check we have at least 2 sequences
    seq_cnt=\$(grep -c '^>' ${query_name}_filtered.fasta || true)
    if [ "\$seq_cnt" -lt 2 ]; then
        echo "WARNING: Only \$seq_cnt sequence(s) after filtering for ${query_name}" >&2
    fi
    """
}

process MAFFT_ALIGN {
    tag "${query_name}"
    publishDir "${params.outdir}/allele_seq", mode: 'copy'

    input:
    tuple val(query_name), path(filtered_fasta)

    output:
    tuple val(query_name), path("${query_name}_aligned.fasta"), optional: true

    script:
    """
    seq_cnt=\$(grep -c '^>' ${filtered_fasta} || true)
    if [ "\$seq_cnt" -ge 2 ]; then
        mafft --auto ${filtered_fasta} > ${query_name}_aligned.fasta
    else
        echo "Skipping MAFFT: only \$seq_cnt sequence(s)" >&2
    fi
    """
}

process WEBLOGO {
    tag "${query_name}"
    publishDir "${params.outdir}/logo", mode: 'copy'

    input:
    tuple val(query_name), path(aligned_fasta)

    output:
    path "${query_name}_logo.${params.weblogo_format}"

    script:
    """
    # Expand sequences by species weight
    awk -v wrap=60 '
        BEGIN { RS=">"; ORS=""; }
        NR > 1 {
            n = split(\$0, L, "\\n")
            header = L[1]
            seq = ""
            for (i = 2; i <= n; i++) seq = seq L[i]

            dup = 1
            if (match(header, /species_number=([0-9]+)/, m)) dup = m[1] + 0
            else if (match(header, /species_percentage=([0-9.]+)/, m))
                dup = int(m[1] + 0.5)

            for (d = 1; d <= dup; d++) {
                printf(">%s_%d\\n", header, d)
                for (p = 1; p <= length(seq); p += wrap)
                    printf("%s\\n", substr(seq, p, wrap))
            }
        }
    ' ${aligned_fasta} > expanded.fasta

    weblogo \
        --fin expanded.fasta \
        --fout ${query_name}_logo.${params.weblogo_format} \
        --format ${params.weblogo_format} \
        --datatype fasta \
        --title "${query_name}" \
        --units ${params.weblogo_units} \
        --color-scheme ${params.weblogo_color_scheme} \
        --size ${params.weblogo_size} \
        --resolution ${params.weblogo_resolution}
    """
}

process ALIGNMENT_FIGURE {
    tag "${query_name}/${ref_type}"
    publishDir "${params.outdir}/alignment_figures/${query_name}/${ref_type}", mode: 'copy'

    input:
    tuple val(query_name), val(ref_type), path(csv_files)

    output:
    tuple val(query_name), path("*.pdf"), optional: true

    script:
    """
    python ${params.code_dir}/alignment_figure.py \
        allele_info_no_strain_name.csv \
        ${query_name}_alignment.pdf \
        ${params.species_threshold} \
    || echo "Skipping no-gap figure for ${query_name}" >&2

    python ${params.code_dir}/alignment_figure.py \
        allele_info_no_strain_name.csv \
        ${query_name}_alignment_gap.pdf \
        ${params.species_threshold} \
        --gap \
    || echo "Skipping gap figure for ${query_name}" >&2
    """
}

process ALIGNMENT_FIGURE_REF {
    tag "${query_name}/${ref_type}"
    publishDir "${params.outdir}/alignment_figures/${query_name}/${ref_type}", mode: 'copy'

    input:
    tuple val(query_name), val(ref_type), path(csv_files)

    output:
    tuple val(query_name), path("*.pdf"), optional: true

    script:
    """
    python ${params.code_dir}/alignment_figure.py \
        allele_info_no_strain_name.csv \
        ${query_name}_alignment.pdf \
        ${params.species_threshold} \
    || echo "Skipping no-gap figure for ${query_name}" >&2

    python ${params.code_dir}/alignment_figure.py \
        allele_info_no_strain_name.csv \
        ${query_name}_alignment_gap.pdf \
        ${params.species_threshold} \
        --gap \
    || echo "Skipping gap figure for ${query_name}" >&2
    """
}
