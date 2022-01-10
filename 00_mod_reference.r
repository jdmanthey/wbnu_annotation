library(Biostrings)

# define output directory
outdir <- "wbnu_scaffolds"
dir.create(outdir)

genome <- readDNAStringSet("WBNU_corrected_nameTrim_final_assembly.fasta")

# keep only scaffolds >= 1Mbp (98.95%, 1,021,383,967 bases)
genome2 <- genome[genome@ranges@width >= 1000000]

genome_names <- genome2@ranges@NAMES[genome2@ranges@width >= 1000000]
genome_output_names <- paste0(outdir, "/", genome_names, ".fasta")

# write complete genome
writeXStringSet(genome2, file="wbnu.fasta")

# write files for each scaffold
for(a in 1:length(genome_names)) {
	a_rep <- genome2[a]
	# write output
	writeXStringSet(a_rep, file=genome_output_names[a])
}
