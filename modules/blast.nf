process MAKE_BLAST_DB {
    tag "${genome.baseName}"

    input:
    path genome

    output:
    tuple path(genome), path("${genome.baseName}.*")

    script:
    def prefix = genome.baseName
    """
    makeblastdb -in ${genome} -dbtype nucl -out ${prefix} -parse_seqids
    """
}

process COLLECT_BLAST_DB {
    tag "${genome.baseName}"

    input:
    path genome

    output:
    tuple path(genome), path("${genome}.n*")

    script:
    """
    # DB files already exist alongside the genome; just symlink them into the work dir
    for f in ${genome.toRealPath()}.n*; do
        ln -s "\$f" . 2>/dev/null || true
    done
    """
}

process BLAST_SEARCH {
    tag "${query.baseName}_vs_${genome.baseName}"
    cpus 4

    input:
    tuple path(query), path(genome), path(db_files)

    output:
    tuple val("${query.baseName}"), path("${query.baseName}_vs_${genome.baseName}.out")

    script:
    def blast_cmd = params.blast_type ?: 'blastn'
    """
    # Determine DB prefix: makeblastdb creates <baseName>.n*, pre-built uses <filename>.n*
    if [ -f "${genome.baseName}.nin" ] || [ -f "${genome.baseName}.nsq" ]; then
        DB_PREFIX="${genome.baseName}"
    else
        DB_PREFIX="${genome.name}"
    fi

    ${blast_cmd} -query ${query} \
        -db \${DB_PREFIX} \
        -out ${query.baseName}_vs_${genome.baseName}.out \
        -evalue 1e-5 \
        -num_threads ${task.cpus}
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
        --similarity 0 \
        --min_length 0 \
        --file_name ${query_name}
    """
}
