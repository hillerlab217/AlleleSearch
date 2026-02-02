#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ──────────────────────────────────────────────────────────────
// Surface Chaperones Allele Discovery Pipeline
// ──────────────────────────────────────────────────────────────

include { MAKE_BLAST_DB; COLLECT_BLAST_DB; BLAST_SEARCH; BLAST_SEARCH_BATCHED; PARSE_BLAST } from './modules/blast'
include { FILTER_BLAST }                             from './modules/filter'
include { CSV_TO_FASTA; CODON_CHECK as CODON_CHECK_10; CODON_CHECK as CODON_CHECK_20; EXTEND_SEQUENCES; COLLECT_STANDARD } from './modules/orf'
include { TRANSLATE }                                from './modules/translate'
include { ALLELE_QUERY_REF; ALLELE_SELF_REF }        from './modules/allele'
include { FILTER_ALLELES; MAFFT_ALIGN; WEBLOGO; TCOFFEE_ALIGN } from './modules/quality'

// ──────────────────────────────────────────────────────────────
// Help message
// ──────────────────────────────────────────────────────────────
if (params.help) {
    log.info """
    =========================================
    Surface Chaperones Allele Pipeline
    =========================================
    Usage:
      nextflow run main.nf --query_dir <path> --genome_dir <path>

    Required:
      --query_dir       Directory with query .fasta files
      --genome_dir      Directory with genome .fasta files
      --genome_ref_dir  Directory with full genomes for upstream/downstream extension

    Optional:
      --code_dir        Path to code/ scripts (default: ../surface_chaperones/code)
      --outdir          Output directory (default: results)
      --blast_type      blastn or tblastn (default: blastn)
      --max_evalue      E-value cutoff (default: 0.00001)
      --similarity      Similarity % cutoff (default: 80)
      --length_percentage  Length % cutoff (default: 80)
      --extension_codons   Codons to extend (default: 10)
      --species_threshold  Min species % for filtering (default: 1)
      --make_blastdb    Build BLAST databases (default: false, uses existing DBs)
      --skip_tcoffee    Skip T-Coffee alignment (default: false)

    Step control (--step):
      all             Run full pipeline (default)
      blast_only      BLAST search only (outputs raw .out files)
      filter          Parse + filter BLAST results (requires --blast_raw_dir)
      post_filter     ORF/translate/allele/align steps (requires --blast_filtered_dir)
      visualization   WebLogo + T-Coffee only (requires --aligned_seq_dir)

    Step-specific inputs:
      --blast_raw_dir      Directory with raw BLAST .out files (for --step filter)
      --blast_filtered_dir Directory with filtered BLAST .csv files (for --step post_filter)
      --aligned_seq_dir    Directory with aligned .fasta files (for --step visualization)
    =========================================
    """.stripIndent()
    exit 0
}

// ──────────────────────────────────────────────────────────────
// Workflow
// ──────────────────────────────────────────────────────────────
workflow {

    def run_blast      = params.step == 'all' || params.step == 'blast_only'
    def run_filter     = params.step == 'all' || params.step == 'filter'
    def run_post       = params.step == 'all' || params.step == 'post_filter'
    def run_viz        = params.step == 'all' || params.step == 'visualization'

    // ── BLAST (steps 1–2) ───────────────────────────────────
    if (run_blast) {
        ch_queries = Channel.fromPath("${params.query_dir}/*.fasta")
        ch_genomes = Channel.fromPath("${params.genome_dir}/*.{fasta,fa,fna,ffn}")

        if (params.make_blastdb) {
            MAKE_BLAST_DB(ch_genomes)
            ch_blast_input = ch_queries
                .combine(MAKE_BLAST_DB.out)
            BLAST_SEARCH(ch_blast_input)
            ch_blast_out = BLAST_SEARCH.out
        } else {
            ch_genome_batches = ch_genomes
                .buffer(size: params.blast_batch_size, remainder: true)
                .map { batch -> [batch] }
            ch_blast_input = ch_queries.combine(ch_genome_batches)
            BLAST_SEARCH_BATCHED(ch_blast_input)
            ch_blast_out = BLAST_SEARCH_BATCHED.out.transpose()
        }
    }

    // ── Filter (steps 3–4) ──────────────────────────────────
    if (run_filter) {
        ch_queries_filter = Channel.fromPath("${params.query_dir}/*.fasta")

        if (params.step == 'filter') {
            // Load raw BLAST .out files from blast_raw_dir
            // File names are like QueryName_vs_GenomeName.out
            ch_blast_out = Channel.fromPath("${params.blast_raw_dir}/*.out")
                .map { f ->
                    def query_name = f.baseName.split('_vs_')[0]
                    [query_name, f]
                }
        }

        ch_blast_grouped = ch_blast_out
            .groupTuple()

        PARSE_BLAST(ch_blast_grouped)

        ch_query_map = ch_queries_filter.map { q -> [q.baseName, q] }
        ch_filter_input = ch_blast_out
            .groupTuple()
            .join(ch_query_map)

        FILTER_BLAST(ch_filter_input)
    }

    // ── Post-filter (steps 5–16) ────────────────────────────
    if (run_post) {
        ch_queries_post = Channel.fromPath("${params.query_dir}/*.fasta")
        ch_query_map_post = ch_queries_post.map { q -> [q.baseName, q] }
        ch_genome_ref = Channel.fromPath(params.genome_ref_dir)

        if (params.step == 'post_filter') {
            // Load filtered BLAST CSVs from blast_filtered_dir
            ch_filter_out = Channel.fromPath("${params.blast_filtered_dir}/*.csv")
                .map { f -> [f.baseName, f] }
        } else {
            ch_filter_out = FILTER_BLAST.out
        }

        // Step 5: CSV to FASTA
        CSV_TO_FASTA(ch_filter_out)

        // Step 6: Check start/stop codons (10-codon window)
        ch_codon10_in = CSV_TO_FASTA.out.map { name, fasta -> [name, fasta, 10, "down"] }
        CODON_CHECK_10(ch_codon10_in)

        // Step 7: Extend sequences ±10 codons from genome
        EXTEND_SEQUENCES(CODON_CHECK_10.out.classified, ch_genome_ref.first())

        // Step 8: Re-check with 20-codon window
        ch_codon20_in = EXTEND_SEQUENCES.out.map { name, fasta -> [name, fasta, 20, "up"] }
        CODON_CHECK_20(ch_codon20_in)

        // Step 9: Collect standard sequences from both rounds
        ch_standard = CODON_CHECK_10.out.standard
            .join(CODON_CHECK_20.out.standard)

        COLLECT_STANDARD(ch_standard)

        // Step 10: Translate DNA to protein
        TRANSLATE(COLLECT_STANDARD.out)

        // Step 11: Allele analysis — query reference
        ch_allele_qref = TRANSLATE.out.inframe
            .join(ch_query_map_post)

        ALLELE_QUERY_REF(ch_allele_qref)

        // Step 12: Allele analysis — self reference
        ALLELE_SELF_REF(TRANSLATE.out.inframe)

        // Step 13: Filter alleles by species percentage
        ch_allele_fasta = ALLELE_SELF_REF.out.results
            .map { name, csv, fasta -> [name, fasta] }

        FILTER_ALLELES(ch_allele_fasta)

        // Step 14: MAFFT alignment
        MAFFT_ALIGN(FILTER_ALLELES.out)

        // Step 15: WebLogo
        WEBLOGO(MAFFT_ALIGN.out)

        // Step 16: T-Coffee (optional)
        if (!params.skip_tcoffee) {
            TCOFFEE_ALIGN(FILTER_ALLELES.out)
        }
    }

    // ── Visualization only (steps 15–16) ────────────────────
    if (params.step == 'visualization') {
        // Load pre-aligned FASTA files
        ch_aligned = Channel.fromPath("${params.aligned_seq_dir}/*_aligned.fasta")
            .map { f ->
                def name = f.baseName.replaceAll(/_aligned$/, '')
                [name, f]
            }

        WEBLOGO(ch_aligned)

        if (!params.skip_tcoffee) {
            // For T-Coffee we need unaligned sequences; use the aligned ones as input
            TCOFFEE_ALIGN(ch_aligned)
        }
    }
}
