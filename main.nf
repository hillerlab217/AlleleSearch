#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ──────────────────────────────────────────────────────────────
// Allele Search Pipeline
// ──────────────────────────────────────────────────────────────

include { COLLECT_ALL_BLAST_DBS; BLAST_SEARCH; PARSE_BLAST } from './modules/blast'
include { FILTER_BLAST }                             from './modules/filter'
include { EXTRACT_DNA; CSV_TO_FASTA; CODON_CHECK as CODON_CHECK_10; CODON_CHECK as CODON_CHECK_20; EXTEND_SEQUENCES; COLLECT_STANDARD } from './modules/orf'
include { TRANSLATE }                                from './modules/translate'
include { ALLELE_SELF_REF; ALLELE_CONSENSUS_REF }    from './modules/allele'
include { MAFFT_FOR_CONSENSUS; COMPUTE_CONSENSUS; FILTER_ALLELES; MAFFT_ALIGN; WEBLOGO; ALIGNMENT_FIGURE; ALIGNMENT_FIGURE_REF } from './modules/quality'

// ──────────────────────────────────────────────────────────────
// Resolve step mode
// ──────────────────────────────────────────────────────────────
def resolveStep() {
    if (params.step)   return params.step
    if (params.blast)  return 'blast'
    if (params.filter) return 'filter'
    if (params.orf)    return 'orf'
    if (params.visual) return 'visual'
    return 'full'
}

def run_step = resolveStep()

def valid_steps = ['full', 'blast', 'filter', 'orf', 'visual']
if (!(run_step in valid_steps)) {
    error "ERROR: Invalid --step '${run_step}'. Valid options: ${valid_steps.join(', ')}"
}

// Check that at most one boolean step flag is set
def step_count = [params.blast, params.filter, params.orf, params.visual].count { it }
if (step_count > 1) {
    error "ERROR: Only one step mode can be specified at a time. Use --blast, --filter, --orf, or --visual."
}

log.info "Running in '${run_step}' mode"

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

    Step modes (run a specific stage only):
      --blast           Run BLAST search + parse + filter
      --filter          Re-filter existing BLAST results (requires prior --blast run)
      --orf             Run ORF finding + translation on existing filtered results
      --visual          Run visualization on existing allele results
      --step <name>     Alternative syntax: blast, filter, orf, visual

    Optional:
      --code_dir        Path to code/ scripts (default: ./code)
      --outdir          Output directory (default: ./results)
      --blast_type      blastn or tblastn (default: blastn)
      --max_evalue      E-value cutoff (default: 0.00001)
      --sequence_identity  Sequence identity % cutoff (default: 80)
      --length_percentage  Length % cutoff (default: 80)
      --extension_codons   Codons to extend (default: 10)
      --species_threshold  Min species % for filtering (default: 1)
      --make_blastdb    Build BLAST databases (default: false, uses existing DBs)

    WebLogo:
      --weblogo_format       Output format: eps, png, pdf, svg, jpeg (default: eps)
      --weblogo_units        Unit type: probability, bits, nats (default: probability)
      --weblogo_color_scheme Color scheme: chemistry, hydrophobicity, charge, classic (default: chemistry)
      --weblogo_size         Logo size: small, medium, large (default: medium)
      --weblogo_resolution   Resolution in DPI (default: 96)
    =========================================
    """.stripIndent()
    exit 0
}

// ──────────────────────────────────────────────────────────────
// Parameter validation (context-dependent)
// ──────────────────────────────────────────────────────────────
if (run_step in ['full', 'blast']) {
    if (!params.query_dir)  { error "ERROR: --query_dir is required for '${run_step}' mode" }
    if (!params.genome_dir) { error "ERROR: --genome_dir is required for '${run_step}' mode" }
}

if (run_step == 'filter') {
    if (!params.query_dir) { error "ERROR: --query_dir is required for 'filter' mode (needed for length calculation)" }
    def blast_summary_dir = file("${params.outdir}/blast_summary")
    if (!blast_summary_dir.isDirectory()) {
        error "ERROR: --filter mode requires existing parsed results in ${params.outdir}/blast_summary/. Run --blast first."
    }
}

if (run_step == 'orf') {
    def filtered_dir = file("${params.outdir}/blast_summary_filtered")
    if (!filtered_dir.isDirectory()) {
        error "ERROR: --orf mode requires existing filtered results in ${params.outdir}/blast_summary_filtered/. Run --blast or --filter first."
    }
    if (!params.genome_dir) { error "ERROR: --genome_dir is required for 'orf' mode (needed for sequence extension)" }
}

if (run_step == 'visual') {
    def allele_filt_dir = file("${params.outdir}/allele_seq_filtered")
    if (!allele_filt_dir.isDirectory()) {
        error "ERROR: --visual mode requires existing allele results in ${params.outdir}/allele_seq_filtered/. Run the full pipeline first."
    }
    def allele_info_dir = file("${params.outdir}/allele_info")
    if (!allele_info_dir.isDirectory()) {
        error "ERROR: --visual mode requires existing allele info in ${params.outdir}/allele_info/. Run the full pipeline first."
    }
}

// ──────────────────────────────────────────────────────────────
// Workflow
// ──────────────────────────────────────────────────────────────
workflow {

    if (run_step == 'blast') {
        // =============================================================
        // BLAST-only mode: DB creation + search + parse + filter
        // =============================================================
        ch_queries = Channel.fromPath("${params.query_dir}/*.{fasta,fa,fna,ffn}")
        ch_genome_dir = Channel.fromPath(params.genome_dir)

        COLLECT_ALL_BLAST_DBS(ch_genome_dir)

        ch_genome_batches = Channel
            .fromPath("${params.genome_dir}/*.{fasta,fa,fna,ffn}")
            .map { it.name }
            .collate(params.blast_batch_size)
            .map { batch -> batch.join(',') }

        ch_blast_input = ch_queries
            .combine(ch_genome_batches)
            .combine(COLLECT_ALL_BLAST_DBS.out)
        BLAST_SEARCH(ch_blast_input)

        ch_grouped_blast = BLAST_SEARCH.out.groupTuple().map { name, files -> [name, files.flatten()] }
        PARSE_BLAST(ch_grouped_blast)

        ch_query_map = ch_queries.map { q -> [q.baseName, q] }
        ch_filter_input = PARSE_BLAST.out.join(ch_query_map)
        FILTER_BLAST(ch_filter_input)

    } else if (run_step == 'filter') {
        // =============================================================
        // Filter-only mode: re-filter existing parsed BLAST CSVs
        // =============================================================
        ch_queries = Channel.fromPath("${params.query_dir}/*.{fasta,fa,fna,ffn}")
        ch_query_map = ch_queries.map { q -> [q.baseName, q] }

        ch_parsed_csv = Channel
            .fromPath("${params.outdir}/blast_summary/*.csv")
            .map { f -> [f.baseName, f] }

        ch_filter_input = ch_parsed_csv.join(ch_query_map)
        FILTER_BLAST(ch_filter_input)

    } else if (run_step == 'orf') {
        // =============================================================
        // ORF-only mode: ORF finding + translation from existing filtered CSVs
        // =============================================================
        ch_genome_ref = Channel.fromPath(params.genome_dir)

        ch_csv_input = Channel
            .fromPath("${params.outdir}/blast_summary_filtered/*.csv")
            .map { f -> [f.baseName, f] }

        if (params.blast_type == 'tblastn') {
            EXTRACT_DNA(ch_csv_input, ch_genome_ref.first())
            ch_csv_for_fasta = EXTRACT_DNA.out
        } else {
            ch_csv_for_fasta = ch_csv_input
        }

        CSV_TO_FASTA(ch_csv_for_fasta)

        ch_codon10_in = CSV_TO_FASTA.out.map { name, fasta -> [name, fasta, 10, "down"] }
        CODON_CHECK_10(ch_codon10_in)

        EXTEND_SEQUENCES(CODON_CHECK_10.out.classified, ch_genome_ref.first())

        ch_codon20_in = EXTEND_SEQUENCES.out.map { name, fasta -> [name, fasta, 20, "up"] }
        CODON_CHECK_20(ch_codon20_in)

        ch_standard = CODON_CHECK_10.out.standard
            .join(CODON_CHECK_20.out.standard)
        COLLECT_STANDARD(ch_standard)

        TRANSLATE(COLLECT_STANDARD.out)

    } else if (run_step == 'visual') {
        // =============================================================
        // Visual-only mode: MAFFT + WebLogo + alignment figures
        // =============================================================

        // MAFFT and WebLogo from filtered allele FASTAs
        ch_filtered_alleles = Channel
            .fromPath("${params.outdir}/allele_seq_filtered/*_filtered.fasta")
            .map { f ->
                def name = f.baseName.replaceAll(/_filtered$/, '')
                [name, f]
            }

        MAFFT_ALIGN(ch_filtered_alleles)
        WEBLOGO(MAFFT_ALIGN.out)

        // Alignment figures from allele info CSVs
        ch_figure_self = Channel
            .fromPath("${params.outdir}/allele_info/*/self_reference/*.csv")
            .map { f ->
                def query_name = f.parent.parent.name
                [query_name, f]
            }
            .groupTuple()
            .map { name, csvs -> [name, "self_reference", csvs] }

        ALIGNMENT_FIGURE(ch_figure_self)

        ch_consensus_csvs = Channel
            .fromPath("${params.outdir}/allele_info/*/consensus_reference/*.csv")
            .map { f ->
                def query_name = f.parent.parent.name
                [query_name, f]
            }
            .groupTuple()
            .map { name, csvs -> [name, "consensus_reference", csvs] }

        ch_consensus_fasta = Channel
            .fromPath("${params.outdir}/allele_info/*/consensus_reference/*_consensus.fasta")
            .map { f -> [f.parent.parent.name, f] }

        ALIGNMENT_FIGURE_REF(ch_consensus_csvs.join(ch_consensus_fasta)
            .map { name, csvs, fasta -> [name, "consensus_reference", csvs, fasta] })

    } else {
        // =============================================================
        // Full pipeline (default)
        // =============================================================
        ch_queries = Channel.fromPath("${params.query_dir}/*.{fasta,fa,fna,ffn}")
        ch_genome_ref = Channel.fromPath(params.genome_dir)

        // Step 1: Ensure all BLAST databases exist (fast directory-level check)
        COLLECT_ALL_BLAST_DBS(ch_genome_ref)

        // Step 2: Run BLAST - batched (100 genomes per task)
        ch_genome_batches = Channel
            .fromPath("${params.genome_dir}/*.{fasta,fa,fna,ffn}")
            .map { it.name }
            .collate(params.blast_batch_size)
            .map { batch -> batch.join(',') }

        ch_blast_input = ch_queries
            .combine(ch_genome_batches)
            .combine(COLLECT_ALL_BLAST_DBS.out)
        BLAST_SEARCH(ch_blast_input)

        // Step 3: Parse BLAST results (collect all batches per query)
        ch_grouped_blast = BLAST_SEARCH.out.groupTuple().map { name, files -> [name, files.flatten()] }
        PARSE_BLAST(ch_grouped_blast)

        // Step 4: Filter BLAST results with strict cutoffs
        ch_query_map = ch_queries.map { q -> [q.baseName, q] }
        ch_filter_input = PARSE_BLAST.out.join(ch_query_map)
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

        // Step 7: Check start/stop codons (10-codon window)
        ch_codon10_in = CSV_TO_FASTA.out.map { name, fasta -> [name, fasta, 10, "down"] }
        CODON_CHECK_10(ch_codon10_in)

        // Step 8: Extend sequences from genome
        EXTEND_SEQUENCES(CODON_CHECK_10.out.classified, ch_genome_ref.first())

        // Step 9: Re-check with 20-codon window
        ch_codon20_in = EXTEND_SEQUENCES.out.map { name, fasta -> [name, fasta, 20, "up"] }
        CODON_CHECK_20(ch_codon20_in)

        // Step 10: Collect standard sequences from both rounds
        ch_standard = CODON_CHECK_10.out.standard
            .join(CODON_CHECK_20.out.standard)
        COLLECT_STANDARD(ch_standard)

        // Step 11: Translate DNA to protein
        TRANSLATE(COLLECT_STANDARD.out)

        // Step 12: Allele analysis - self reference
        ALLELE_SELF_REF(TRANSLATE.out.inframe)

        // Step 13: MAFFT on all alleles → consensus sequence
        MAFFT_FOR_CONSENSUS(ALLELE_SELF_REF.out.results.map { name, csv, fasta -> [name, fasta] })
        COMPUTE_CONSENSUS(MAFFT_FOR_CONSENSUS.out)

        // Step 14: Allele analysis - consensus reference
        ch_allele_cref = TRANSLATE.out.inframe
            .join(COMPUTE_CONSENSUS.out)
        ALLELE_CONSENSUS_REF(ch_allele_cref)

        // Step 15: Filter alleles by species percentage
        FILTER_ALLELES(ALLELE_SELF_REF.out.results.map { name, csv, fasta -> [name, fasta] })

        // Step 16: Alignment variant figures (both reference types)
        ch_figure_self = ALLELE_SELF_REF.out.results
            .map { name, csvs, fasta -> [name, "self_reference", csvs] }
        ALIGNMENT_FIGURE(ch_figure_self)

        ch_figure_consensus = ALLELE_CONSENSUS_REF.out.results
            .map { name, csvs, allele_fasta, cons_fasta -> [name, "consensus_reference", csvs, cons_fasta] }
        ALIGNMENT_FIGURE_REF(ch_figure_consensus)

        // Step 17: MAFFT alignment (filtered alleles for WebLogo)
        MAFFT_ALIGN(FILTER_ALLELES.out)

        // Step 18: WebLogo
        WEBLOGO(MAFFT_ALIGN.out)
    }
}
