

1	main.nf	--visual mode produces no figures on old results    fix: when run on old results give warning to let the user run the whole thing again.
2	Correctness	protein_allele_variant.py, compute_consensus.py	Variant CSV references different consensus than exported FASTA  fix: do not recomput and add consensus sequence as the first sequence in alignment tabel in allele_info
3	alignment_figure.py, compute_consensus.py	Silent figure skip on empty consensus FASTA when there is no consensus sequence, make an figure with self-reference and tell user that in nexflow run reults
4	protein_allele_variant.py	* deletion notation ambiguous with stop codon   fix: change * to - for deletion and insersion
5	start_stop_codon_dna_seq.py	Short-sequence internal stop check regression   fix: give user warning sequence is too short that this will happen
6	ignore
7	Low	alignment_figure.py, protein_allele_variant.py	Duplicate consensus logic   fix: merge them
8	ignore
