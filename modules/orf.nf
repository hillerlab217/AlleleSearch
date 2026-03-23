process EXTRACT_DNA {
    tag "${query_name}"

    input:
    tuple val(query_name), path(csv_file)
    path genome_dir

    output:
    tuple val(query_name), path("${query_name}_dna.csv")

    script:
    """
    python ${params.code_dir}/extract_dna_from_blast_csv.py \
        ${csv_file} ${query_name}_dna.csv ${genome_dir}
    """
}

process CSV_TO_FASTA {
    tag "${query_name}"

    input:
    tuple val(query_name), path(csv_file)

    output:
    tuple val(query_name), path("${query_name}.fasta")

    script:
    """
    python ${params.code_dir}/csv_to_fasta.py ${csv_file} ${query_name}.fasta
    """
}

process CODON_CHECK {
    tag "${query_name}_${codon_range}codon"

    input:
    tuple val(query_name), path(fasta_file), val(codon_range), val(direction)

    output:
    tuple val(query_name), val(codon_range), path("codon_out/*"), emit: classified
    tuple val(query_name), path("standard_${codon_range}.fasta"), optional: true, emit: standard

    script:
    """
    mkdir -p codon_out
    python ${params.code_dir}/start_stop_codon_dna_seq.py \
        ${fasta_file} \
        codon_out \
        --codon_range ${codon_range} \
        --start_search_direction ${direction}
    # Ensure standard.fasta exists even if empty
    touch codon_out/standard.fasta
    cp codon_out/standard.fasta standard_${codon_range}.fasta
    """
}

process EXTEND_SEQUENCES {
    tag "${query_name}"

    input:
    tuple val(query_name), val(codon_range), path(classified_files)
    path genome_dir

    output:
    tuple val(query_name), path("extended.fasta")

    script:
    def extension_bp = params.extension_codons * 3
    """
    > extended.fasta

    for file in ${classified_files}; do
        filename=\$(basename "\$file")
        run_mode=""
        case "\$filename" in
            no_stop.fasta)          run_mode="down_stream" ;;
            no_start.fasta)         run_mode="up_stream" ;;
            no_stop_no_start.fasta) run_mode="up_down_stream" ;;
            *)                      continue ;;
        esac

        python ${params.code_dir}/up_down_stream.py \
            --mode "\$run_mode" \
            --seq_file "\$file" \
            --genome_dir ${genome_dir} \
            --length ${extension_bp} \
            --output_fasta temp_ext.fasta || true

        if [ -f temp_ext.fasta ]; then
            cat temp_ext.fasta >> extended.fasta
            rm temp_ext.fasta
        fi
    done
    """
}

process COLLECT_STANDARD {
    tag "${query_name}"
    publishDir "${params.outdir}/standard_dna_seq", mode: 'copy'

    input:
    tuple val(query_name), path(standard_10), path(standard_20)

    output:
    tuple val(query_name), path("${query_name}.fasta"), optional: true

    script:
    """
    cat ${standard_10} ${standard_20} > ${query_name}.fasta
    seq_cnt=\$(grep -c '^>' ${query_name}.fasta || echo 0)
    if [ "\$seq_cnt" -eq 0 ]; then
        rm ${query_name}.fasta
        echo "WARNING: [${query_name}] No valid sequences after ORF extension — skipping downstream steps"
    fi
    """
}
