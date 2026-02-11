#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ──────────────────────────────────────────────────────────────
// Allele Search Pipeline
// ──────────────────────────────────────────────────────────────

include { COLLECT_BLAST_DB; BLAST_SEARCH; PARSE_BLAST } from './modules/blast'
include { FILTER_BLAST }                             from './modules/filter'
include { EXTRACT_DNA; CSV_TO_FASTA; CODON_CHECK as CODON_CHECK_10; CODON_CHECK as CODON_CHECK_20; EXTEND_SEQUENCES; COLLECT_STANDARD } from './modules/orf'
include { TRANSLATE }                                from './modules/translate'
include { ALLELE_QUERY_REF; ALLELE_SELF_REF }        from './modules/allele'
include { FILTER_ALLELES; MAFFT_ALIGN; WEBLOGO; ALIGNMENT_FIGURE; TCOFFEE_ALIGN } from './modules/quality'

// ──────────────────────────────────────────────────────────────
// Help message
// ──────────────────────────────────────────────────────────────
if (params.help) {
    log.info """
    =========================================
    Allele Search Pipeline
    =========================================
    Usage:
      nextflow run main.nf --query_dir <path> --genome_dir <path>

    Required:
      --query_dir       Directory with query .fasta files
      --genome_dir      Directory with genome .fasta files (also used for sequence extension)

    Optional:
      --code_dir        Path to code/ scripts (default: ./code)
      --outdir          Output directory (default: ./results)
      --blast_type      blastn or tblastn (default: blastn)
      --max_evalue      E-value cutoff (default: 0.00001)
      --similarity      Similarity % cutoff (default: 80)
      --length_percentage  Length % cutoff (default: 80)
      --extension_codons   Codons to extend (default: 10)
      --species_threshold  Min species % for filtering (default: 1)
      --make_blastdb    Build BLAST databases (default: false, uses existing DBs)
      --skip_tcoffee    Skip T-Coffee alignment (default: false)
    =========================================
    """.stripIndent()
    exit 0
}

// ──────────────────────────────────────────────────────────────
// Parameter validation
// ──────────────────────────────────────────────────────────────
if (!params.query_dir)      { error "ERROR: --query_dir is required" }
if (!params.genome_dir)     { error "ERROR: --genome_dir is required" }

// ──────────────────────────────────────────────────────────────
// Input channels
// ──────────────────────────────────────────────────────────────
ch_queries = Channel.fromPath("${params.query_dir}/*.{fasta,fa,fna,ffn}")
ch_genomes = Channel.fromPath("${params.genome_dir}/*.{fasta,fa,fna,ffn}")
ch_genome_ref = Channel.fromPath(params.genome_dir)

// ──────────────────────────────────────────────────────────────
// Workflow
// ──────────────────────────────────────────────────────────────
workflow {

    // Step 1: Collect pre-built BLAST databases or build them on-the-fly
    COLLECT_BLAST_DB(ch_genomes)
    ch_db = COLLECT_BLAST_DB.out

    // Step 2: Run BLAST — each query × each genome
    ch_blast_input = ch_queries
        .combine(ch_db)  // [query, genome, db_files]

    BLAST_SEARCH(ch_blast_input)

    // Step 3: Parse BLAST results (group .out files by query name)
    ch_blast_grouped = BLAST_SEARCH.out
        .groupTuple()  // [query_name, [out_files]]

    PARSE_BLAST(ch_blast_grouped)

    // Step 4: Filter BLAST results with strict cutoffs
    ch_query_map = ch_queries.map { q -> [q.baseName, q] }
    ch_filter_input = BLAST_SEARCH.out
        .groupTuple()
        .join(ch_query_map)  // [query_name, [blast_files], query_fasta]

    FILTER_BLAST(ch_filter_input)

    // Step 5: For tblastn, extract DNA from genomes using subject coordinates
    if (params.blast_type == 'tblastn') {
        EXTRACT_DNA(FILTER_BLAST.out, ch_genome_ref.first())
        ch_csv_for_fasta = EXTRACT_DNA.out
    } else {
        ch_csv_for_fasta = FILTER_BLAST.out
    }

    // Step 6: CSV to FASTA
    CSV_TO_FASTA(ch_csv_for_fasta)

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
        .join(CODON_CHECK_20.out.standard)  // [query_name, std_10, std_20]

    COLLECT_STANDARD(ch_standard)

    // Step 10: Translate DNA to protein
    TRANSLATE(COLLECT_STANDARD.out)

    // Step 11: Allele analysis — query reference
    ch_allele_qref = TRANSLATE.out.inframe
        .join(ch_query_map)  // [query_name, protein_fasta, query_fasta]

    ALLELE_QUERY_REF(ch_allele_qref)

    // Step 12: Allele analysis — self reference
    ALLELE_SELF_REF(TRANSLATE.out.inframe)

    // Step 13: Filter alleles by species percentage
    ch_allele_fasta = ALLELE_SELF_REF.out.results
        .map { name, csv, fasta -> [name, fasta] }

    FILTER_ALLELES(ch_allele_fasta)

    // Step 14: Alignment variant figures (both reference types)
    ch_figure_self = ALLELE_SELF_REF.out.results
        .map { name, csvs, fasta -> [name, "self_reference", csvs] }
    ch_figure_query = ALLELE_QUERY_REF.out
        .map { name, csvs, fasta -> [name, "query_reference", csvs] }

    ALIGNMENT_FIGURE(ch_figure_self.mix(ch_figure_query))

    // Step 15: MAFFT alignment
    MAFFT_ALIGN(FILTER_ALLELES.out)

    // Step 16: WebLogo
    WEBLOGO(MAFFT_ALIGN.out)

    // Step 17: T-Coffee (optional)
    if (!params.skip_tcoffee) {
        TCOFFEE_ALIGN(FILTER_ALLELES.out)
    }
}
