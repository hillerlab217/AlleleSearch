process COLLECT_ALL_BLAST_DBS {
    tag "all_genomes"

    input:
    path genome_dir

    output:
    path genome_dir

    script:
    """
    if [ "${params.make_blastdb}" = "true" ]; then
        # Find genomes missing BLAST DBs (directory-level diff, not per-file ls)
        MISSING=\$(comm -23 \
            <(ls ${genome_dir}/*.fasta ${genome_dir}/*.fa ${genome_dir}/*.fna ${genome_dir}/*.ffn 2>/dev/null | sort) \
            <(ls ${genome_dir}/*.nin 2>/dev/null | sed 's/\\.nin\$//' | sort))

        if [ -n "\$MISSING" ]; then
            echo "Building BLAST DBs for \$(echo "\$MISSING" | wc -l) genomes..."
            for genome in \$MISSING; do
                makeblastdb -in "\$genome" -dbtype nucl -out "\$genome" -parse_seqids
            done
        else
            echo "All BLAST DBs already exist, skipping."
        fi
    else
        echo "make_blastdb=false: skipping DB creation."
    fi
    """
}

process BLAST_SEARCH {
    tag "${query.baseName}_batch"
    cpus 2
    publishDir "${params.outdir}/blast_raw", mode: 'copy'

    input:
    tuple path(query), val(genome_names_csv), path(genome_dir)

    output:
    tuple val("${query.baseName}"), path("*.out")

    script:
    def blast_cmd = params.blast_type ?: 'blastn'
    def shared_args = [
        params.blast_word_size != null ? "-word_size ${params.blast_word_size}" : '',
        params.blast_gapopen   != null ? "-gapopen ${params.blast_gapopen}"     : '',
        params.blast_gapextend != null ? "-gapextend ${params.blast_gapextend}" : '',
    ].findAll { it }.join(' ')
    def tblastn_args = (blast_cmd == 'tblastn') ? [
        params.tblastn_matrix     ? "-matrix ${params.tblastn_matrix}"                        : '',
        params.tblastn_threshold  != null ? "-threshold ${params.tblastn_threshold}"          : '',
        params.tblastn_comp_based != null ? "-comp_based_stats ${params.tblastn_comp_based}"  : '',
        params.tblastn_seg        ? "-seg ${params.tblastn_seg}"                              : '',
    ].findAll { it }.join(' ') : ''
    """
    for gname in \$(echo '${genome_names_csv}' | tr ',' ' '); do
        genome="${genome_dir}/\${gname}"
        [ ! -e "\$genome" ] && continue
        base=\$(echo "\$gname" | sed 's/\\.[^.]*\$//')

        if [ -f "\${genome}.nin" ] || [ -f "\${genome}.nsq" ]; then
            DB_PREFIX="\${genome}"
            USE_DB=true
        elif [ -f "\${genome%.*}.nin" ] || [ -f "\${genome%.*}.nsq" ]; then
            DB_PREFIX="\${genome%.*}"
            USE_DB=true
        else
            USE_DB=false
        fi

        if [ "\$USE_DB" = "true" ]; then
            ${blast_cmd} -query ${query} \
                -db \${DB_PREFIX} \
                -out ${query.baseName}_vs_\${base}.out \
                -evalue ${params.max_evalue} \
                ${shared_args} ${tblastn_args} \
                -num_threads ${task.cpus}
        else
            ${blast_cmd} -query ${query} \
                -subject \${genome} \
                -out ${query.baseName}_vs_\${base}.out \
                -evalue ${params.max_evalue} \
                ${shared_args} ${tblastn_args}
        fi
    done
    """

    stub:
    """
    touch ${query.baseName}_vs_stub.out
    """
}

process PARSE_BLAST {
    tag "${query_name}"
    publishDir "${params.outdir}/blast_summary", mode: 'copy'

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
        --blast_type ${params.blast_type ?: 'blastn'} \
        --max_Eval 10 \
        --sequence_identity 0 \
        --min_length 0 \
        --file_name ${query_name}
    """

    stub:
    """
    echo "ID,Contig,E-value,Alignment Length,Sequence Identity,Query Start,Subject Start,Subject End,Sequence" > ${query_name}.csv
    """
}
