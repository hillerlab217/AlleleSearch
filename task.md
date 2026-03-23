Phase 1 - Done
1. Modify the visualization step to completely remove the T-Coffee alignment track. Ensure the remaining tracks adjust their layout accordingly.
2. Identify all adjustable parameters for the Sequence Logo generation. Expose these as params in Nextflow so they can be adjusted without changing the source code.
3. add --step parameter for partical pipline run/specific functional run. 
Important[When these specific steps are asked to run, if output already exist rewrite it. if -resume is not specified, reults should be over writen.]
--blast: Run BLAST only. Default behavior: No filtering or parsing. If a filter is specified in args, apply it; otherwise, return raw results.
--filter: Skip BLAST (assume output exists, if not output error). Run the parsing and filtering logic on existing BLAST results. Requires filter parameters (sequence_identity, length, e-value) to be provided.
--visual: Run the visualization step only using existing intermediate files maff file (/allele_info).
--orf: Run the ORF finding step only.(assuming needed previous steps exits)
4. do tests verify that each of the --step modes works correctly and that the full pipeline still runs end-to-end.
5. update the conda yml file as needed, update the readme.md for these changes, make other documentation for best practice


Phase 2 — Fix soon (robustness)
5. Extract FASTA header parsing into a shared utility
6. Move filter.nf Python logic to code/filter_blast.py
7. Make orf.nf handle the empty-extension case properly
8. Strengthen contig matching in extract_dna_from_blast_csv.py

Phase 3 — Nice to have (style/maintainability)
9. Fix typo allele_sufix → allele_suffix
10. Remove duplicate imports
11. Remove import glob dead import
12. Decide on strain_sequence_analysis.py — integrate or document

Want me to start fixing any of these?