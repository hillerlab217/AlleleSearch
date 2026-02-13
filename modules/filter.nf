process FILTER_BLAST {
    tag "${query_name}"
    publishDir "${params.outdir}/blast_summary_filtered", mode: 'copy'

    input:
    tuple val(query_name), path('input.csv'), path(query_fasta)

    output:
    tuple val(query_name), path("${query_name}.csv")

    script:
    """
    # Calculate min_length from query sequence length
    query_len=\$(awk 'BEGIN{seq=""} !/^>/ {seq=seq\$0} END{print length(seq)}' ${query_fasta})
    min_len=\$((query_len * ${params.length_percentage} / 100))

    python3 -c "
import csv, sys

max_ev  = float(sys.argv[1])
sim     = float(sys.argv[2])
min_l   = int(sys.argv[3])
infile  = sys.argv[4]
outfile = sys.argv[5]

with open(infile) as fin, open(outfile, 'w', newline='') as fout:
    reader = csv.reader(fin)
    writer = csv.writer(fout)
    header = next(reader)
    writer.writerow(header)
    ei = header.index('E-value')
    li = header.index('Alignment Length')
    pi = header.index('Sequence Identity') if 'Sequence Identity' in header else header.index('Alignment Percentage')
    for row in reader:
        if float(row[ei]) < max_ev and float(row[pi]) >= sim and int(row[li]) >= min_l:
            writer.writerow(row)
" ${params.max_evalue} ${params.sequence_identity} \$min_len input.csv ${query_name}.csv
    """
}
