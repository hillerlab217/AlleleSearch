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
include { BLAST_SEARCH_DNA; PARSE_BLAST_DNA; FILTER_BLAST_DNA; MAFFT_DNA; DNA_CONSERVATION_FIGURE } from './modules/dna_variation'
include { EXTRACT_NEIGHBOURS }                               from './modules/neighbour'

// ──────────────────────────────────────────────────────────────
// Resolve step mode
// ──────────────────────────────────────────────────────────────
def resolveStep() {
    if (params.step)          return params.step
    if (params.blast)         return 'blast'
    if (params.filter)        return 'filter'
    if (params.orf)           return 'orf'
    if (params.visual)        return 'visual'
    if (params.dna_variation) return 'dna_variation'
    if (params.neighbour)     return 'neighbour'
    return 'full'
}

def run_step = resolveStep()

def valid_steps = ['full', 'blast', 'filter', 'orf', 'visual', 'dna_variation', 'neighbour']
if (!(run_step in valid_steps)) {
    error "ERROR: Invalid --step '${run_step}'. Valid options: ${valid_steps.join(', ')}"
}

// Check that at most one boolean step flag is set
def step_count = [params.blast, params.filter, params.orf, params.visual, params.dna_variation, params.neighbour].count { it }
if (step_count > 1) {
    error "ERROR: Only one step mode can be specified at a time. Use --blast, --filter, --orf, --visual, or --neighbour."
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
      --query_dir       Directory with query .fasta files, or a single query file
      --genome_dir      Directory with genome .fasta files (also used for sequence extension)

    Step modes (run a specific stage only):
      --blast           Run BLAST search + parse + filter
      --filter          Re-filter existing BLAST results (requires prior --blast run)
      --orf             Run ORF finding + translation on existing filtered results
      --visual          Run visualization on existing allele results
      --dna_variation   Run DNA conservation figure (always uses blastn).
                        With --query_dir and --genome_dir: runs BLAST automatically.
                        Without them: reuses existing dna_variation/blast_summary_filtered/ CSVs.
      --step <name>     Alternative syntax: blast, filter, orf, visual, dna_variation

    Optional:
      --code_dir        Path to code/ scripts (default: ./code)
      --outdir          Output directory (default: ./results)
      --blast_type      blastn or tblastn (default: blastn)
      --max_evalue      E-value cutoff (default: 0.00001)
      --sequence_identity  Sequence identity % cutoff (default: 80)
      --length_percentage  Length % cutoff (default: 80)
      --extension_codons   Codons to extend (default: 10)
      --start_codons       Comma-separated valid start codons (default: ATG,GTG,TTG)
      --stop_codons        Comma-separated valid stop codons  (default: TAA,TGA,TAG)
                           Note: non-ATG starts (e.g. GTG, TTG) use a length-based
                           heuristic in translation since protein M-detection won't apply.
      --species_threshold  Min species % for filtering (default: 1)
      --make_blastdb    Build missing BLAST databases before searching (default: false).
                        When false and no DB exists, uses -subject <genome_file> directly.

    BLAST tuning — shared (blastn + tblastn):
      --blast_word_size N       Word size for seeding (default: BLAST default)
      --blast_gapopen N         Gap open penalty (default: BLAST default)
      --blast_gapextend N       Gap extend penalty (default: BLAST default)

    BLAST tuning — tblastn only (ignored when --blast_type blastn):
      --tblastn_matrix NAME     Substitution matrix, e.g. PAM30, BLOSUM45 (default: BLOSUM62)
      --tblastn_threshold N     Minimum word score threshold (default: BLAST default)
      --tblastn_comp_based N    Composition-based stats: 0=off, 1, 2 (default: BLAST default)
      --tblastn_seg VALUE       SEG low-complexity filter: yes or no (default: BLAST default)

      --self_reference_figures  Also generate self-reference alignment figures (default: false)
      --dna_variation_figure_format  Output format for conservation figure: pdf, png, svg (default: pdf)

    WebLogo:
      --weblogo_format       Output format: eps, png, pdf, svg, jpeg (default: eps)
      --weblogo_units        Unit type: probability, bits, nats (default: probability)
      --weblogo_color_scheme Color scheme: chemistry, hydrophobicity, charge, classic (default: chemistry)
      --weblogo_size         Logo size: small, medium, large (default: medium)
      --weblogo_resolution   Resolution in DPI (default: 96)

    Neighbour extraction (runs after visualization when --upstream or --downstream is given):
      --neighbour           Standalone step mode: extract neighbours from existing standard_dna_seq/
      --upstream N          bp upstream of gene to extract (default: 0, 0 = skip)
      --downstream N        bp downstream of gene to extract (default: 0, 0 = skip)
      --flank_only          Exclude the gene itself; output flanking regions only (default: include gene)
      --neighbour_merge     Merge upstream+gene+downstream into one file (requires include_gene; default: split)
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

// Neighbour extraction param validation (applies to both full pipeline and --neighbour step mode)
if (params.flank_only && params.neighbour_merge) {
    error "ERROR: --neighbour_merge requires the gene to be included. --flank_only cannot be combined with --neighbour_merge. Use --split (default) instead."
}

if (run_step == 'neighbour') {
    if (params.upstream == 0 && params.downstream == 0) {
        error "ERROR: --neighbour requires at least --upstream N or --downstream N."
    }
    if (!params.genome_dir) {
        error "ERROR: --neighbour requires --genome_dir."
    }
    def std_dir = file("${params.outdir}/standard_dna_seq")
    if (!std_dir.isDirectory()) {
        error "ERROR: --neighbour mode requires existing sequences in ${params.outdir}/standard_dna_seq/. Run the full pipeline first."
    }
}

if (run_step == 'dna_variation') {
    def has_inputs = params.query_dir && params.genome_dir
    if (!has_inputs) {
        // Reuse-only mode: prior filtered results must exist
        def filtered_dir = file("${params.outdir}/dna_variation/blast_summary_filtered")
        if (!filtered_dir.isDirectory()) {
            error "ERROR: --dna_variation without --query_dir/--genome_dir requires existing filtered results in ${params.outdir}/dna_variation/blast_summary_filtered/. Provide --query_dir and --genome_dir to run BLAST, or run --blast first."
        }
    }
    // blast_type is ignored — dna_variation always uses blastn internally
    if (params.blast_type == 'tblastn') {
        log.warn "WARNING: --blast_type tblastn is ignored for --dna_variation; blastn will be used."
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
        ch_queries = file(params.query_dir).isFile()
            ? Channel.fromPath(params.query_dir)
            : Channel.fromPath("${params.query_dir}/*.{fasta,fa,fna,ffn}")
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
        ch_queries = file(params.query_dir).isFile()
            ? Channel.fromPath(params.query_dir)
            : Channel.fromPath("${params.query_dir}/*.{fasta,fa,fna,ffn}")
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

        ch_codon10_in = CSV_TO_FASTA.out.map { name, fasta -> [name, fasta, params.extension_codons, "down", "r1", params.start_codons, params.stop_codons] }
        CODON_CHECK_10(ch_codon10_in)

        EXTEND_SEQUENCES(CODON_CHECK_10.out.classified, ch_genome_ref.first())

        ch_codon20_in = EXTEND_SEQUENCES.out.map { name, fasta -> [name, fasta, params.extension_codons, "up", "r2", params.start_codons, params.stop_codons] }
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
        if (params.self_reference_figures) {
            ch_figure_self = Channel
                .fromPath("${params.outdir}/allele_info/*/self_reference/*.csv")
                .map { f ->
                    def query_name = f.parent.parent.name
                    [query_name, f]
                }
                .groupTuple()
                .map { name, csvs -> [name, "self_reference", csvs] }
            ALIGNMENT_FIGURE(ch_figure_self)
        }

        ch_consensus_csvs = Channel
            .fromPath("${params.outdir}/allele_info/*/consensus_reference/*.csv")
            .map { f ->
                def query_name = f.parent.parent.name
                [query_name, f]
            }
            .groupTuple()
            .map { name, csvs -> [name, "consensus_reference", csvs] }

        ALIGNMENT_FIGURE_REF(ch_consensus_csvs)

    } else if (run_step == 'dna_variation') {
        // =============================================================
        // DNA variation mode: BLAST (blastn) + filter + MAFFT + conservation figure
        // If --query_dir and --genome_dir are provided: run BLAST automatically.
        // Otherwise: reuse existing dna_variation/blast_summary_filtered/ CSVs.
        // =============================================================

        if (params.query_dir && params.genome_dir) {
            // ── Run BLAST (always blastn) ──────────────────────────────
            ch_queries_dv = file(params.query_dir).isFile()
                ? Channel.fromPath(params.query_dir)
                : Channel.fromPath("${params.query_dir}/*.{fasta,fa,fna,ffn}")
            ch_genome_dir_dv = Channel.fromPath(params.genome_dir)

            COLLECT_ALL_BLAST_DBS(ch_genome_dir_dv)

            ch_genome_batches_dv = Channel
                .fromPath("${params.genome_dir}/*.{fasta,fa,fna,ffn}")
                .map { it.name }
                .collate(params.blast_batch_size)
                .map { batch -> batch.join(',') }

            ch_blast_input_dv = ch_queries_dv
                .combine(ch_genome_batches_dv)
                .combine(COLLECT_ALL_BLAST_DBS.out)
            BLAST_SEARCH_DNA(ch_blast_input_dv)

            ch_grouped_dv = BLAST_SEARCH_DNA.out.groupTuple().map { name, files -> [name, files.flatten()] }
            PARSE_BLAST_DNA(ch_grouped_dv)

            ch_query_map_dv = ch_queries_dv.map { q -> [q.baseName, q] }
            ch_filter_input_dv = PARSE_BLAST_DNA.out.join(ch_query_map_dv)
            FILTER_BLAST_DNA(ch_filter_input_dv)

            ch_filtered_dv = FILTER_BLAST_DNA.out

        } else {
            // ── Reuse existing filtered results ────────────────────────
            ch_filtered_dv = Channel
                .fromPath("${params.outdir}/dna_variation/blast_summary_filtered/*.csv")
                .map { f -> [f.baseName, f] }
        }

        CSV_TO_FASTA(ch_filtered_dv)
        MAFFT_DNA(CSV_TO_FASTA.out)
        DNA_CONSERVATION_FIGURE(MAFFT_DNA.out)

    } else if (run_step == 'neighbour') {
        // =============================================================
        // Neighbour-only mode: extract flanking regions from existing standard_dna_seq/
        // =============================================================
        log.info """
        Neighbour extraction:
          Upstream:   ${params.upstream   > 0 ? params.upstream   + ' bp' : 'not requested (0)'}
          Downstream: ${params.downstream > 0 ? params.downstream + ' bp' : 'not requested (0)'}
          Gene:       ${params.flank_only ? 'excluded (--flank_only)' : 'included (default)'}
          Layout:     ${params.neighbour_merge ? 'merged \u2192 merged.fasta' : 'split \u2192 upstream.fasta / downstream.fasta'}
        """.stripIndent()

        ch_standard_nb = Channel
            .fromPath("${params.outdir}/standard_dna_seq/*.fasta")
            .map { f -> [f.baseName, f] }
        EXTRACT_NEIGHBOURS(ch_standard_nb, params.genome_dir)

    } else {
        // =============================================================
        // Full pipeline (default)
        // =============================================================
        ch_queries = file(params.query_dir).isFile()
            ? Channel.fromPath(params.query_dir)
            : Channel.fromPath("${params.query_dir}/*.{fasta,fa,fna,ffn}")
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

        // Step 7: Check start/stop codons
        ch_codon10_in = CSV_TO_FASTA.out.map { name, fasta -> [name, fasta, params.extension_codons, "down", "r1", params.start_codons, params.stop_codons] }
        CODON_CHECK_10(ch_codon10_in)

        // Step 8: Extend sequences from genome
        EXTEND_SEQUENCES(CODON_CHECK_10.out.classified, ch_genome_ref.first())

        // Step 9: Re-check with same codon window on extended sequences
        ch_codon20_in = EXTEND_SEQUENCES.out.map { name, fasta -> [name, fasta, params.extension_codons, "up", "r2", params.start_codons, params.stop_codons] }
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
        ch_allele_fasta_nonempty = ALLELE_SELF_REF.out.results
            .map { name, csv, fasta -> [name, fasta] }
            .filter { name, fasta -> fasta.text.contains('>') }
        MAFFT_FOR_CONSENSUS(ch_allele_fasta_nonempty)
        COMPUTE_CONSENSUS(MAFFT_FOR_CONSENSUS.out)

        // Step 14: Allele analysis - consensus reference (uses MAFFT alignment for indel-aware variant calling)
        ch_allele_cref = TRANSLATE.out.inframe
            .join(COMPUTE_CONSENSUS.out)
            .join(MAFFT_FOR_CONSENSUS.out)
        ALLELE_CONSENSUS_REF(ch_allele_cref)

        // Step 15: Filter alleles by species percentage
        FILTER_ALLELES(ALLELE_SELF_REF.out.results.map { name, csv, fasta -> [name, fasta] })

        // Step 16: Alignment variant figures
        if (params.self_reference_figures) {
            ch_figure_self = ALLELE_SELF_REF.out.results
                .map { name, csvs, fasta -> [name, "self_reference", csvs] }
            ALIGNMENT_FIGURE(ch_figure_self)
        }

        ch_figure_consensus = ALLELE_CONSENSUS_REF.out.results
            .map { name, csvs, allele_fasta, cons_fasta -> [name, "consensus_reference", csvs] }
        ALIGNMENT_FIGURE_REF(ch_figure_consensus)

        // Step 17: MAFFT alignment (filtered alleles for WebLogo)
        MAFFT_ALIGN(FILTER_ALLELES.out)

        // Step 18: WebLogo
        WEBLOGO(MAFFT_ALIGN.out)

        // Step 19 (optional): Neighbour region extraction
        if (params.upstream > 0 || params.downstream > 0) {
            log.info """
            Neighbour extraction:
              Upstream:   ${params.upstream   > 0 ? params.upstream   + ' bp' : 'not requested (0)'}
              Downstream: ${params.downstream > 0 ? params.downstream + ' bp' : 'not requested (0)'}
              Gene:       ${params.flank_only ? 'excluded (--flank_only)' : 'included (default)'}
              Layout:     ${params.neighbour_merge ? 'merged \u2192 merged.fasta' : 'split \u2192 upstream.fasta / downstream.fasta'}
            """.stripIndent()

            EXTRACT_NEIGHBOURS(COLLECT_STANDARD.out, params.genome_dir)
        }
    }
}
