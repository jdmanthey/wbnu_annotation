options(scipen=999)

# define window size
window_size <- 100000

# output name
output_name <- "wbnu_CDS_TE_summary_100kbp.txt"

# read in repeatmasker output (with overlaps removed)
x <- read.table("wbnu.fasta.no_dups.out", header=T, stringsAsFactors=F, fill=T)
cds <- read.table("wbnu_round2_cds.gff", header=F, stringsAsFactors=F)

# extract just CR1s and ERVs (most common types in genome)
ERV <- x[grepl("LTR/ERV", x$class.family), ]
CR1 <- x[grepl("LINE/CR1", x$class.family), ]

# read in genome index file
genome_index <- read.table("wbnu.fasta.fai", stringsAsFactors=F)

# determine windows
for(a in 1:nrow(genome_index)) {
	n_windows <- floor(genome_index[a,2] / window_size)
	if(a == 1) {
		window_matrix <- data.frame(scaffold=rep(genome_index[a,1], n_windows), start=as.numeric(seq(from=1, to=n_windows) * window_size - window_size + 1), end=as.numeric(seq(from=1, to=n_windows) * window_size))
	} else {
		a_temp <- data.frame(scaffold=rep(genome_index[a,1], n_windows), start=as.numeric(seq(from=1, to=n_windows) * window_size - window_size + 1), end=as.numeric(seq(from=1, to=n_windows) * window_size))
		window_matrix <- rbind(window_matrix, a_temp)
	}
}

# loop for each scaffold to get lists of TE and CDS locations
TE_pos <- list()
ERV_pos <- list()
CR1_pos <- list()
CDS_pos <- list()
for(a in 1:nrow(genome_index)) {
	print(a)
	# get locations of all TEs per scaffold
	a_temp <- x[x$query == genome_index[a,1], ]
	temp_pos <- list()
	for(b in 1:nrow(a_temp)) {
		temp_pos[[b]] <- seq(from=a_temp$q_start[b], to=a_temp$q_end[b], by=1)
	}
	TE_pos[[a]] <- unlist(temp_pos)
	
	# get locations of all CR1s per scaffold
	a_temp <- CR1[CR1$query == genome_index[a,1], ]
	temp_pos <- list()
	for(b in 1:nrow(a_temp)) {
		temp_pos[[b]] <- seq(from=a_temp$q_start[b], to=a_temp$q_end[b], by=1)
	}
	CR1_pos[[a]] <- unlist(temp_pos)

	# get locations of all ERVs per scaffold
	a_temp <- ERV[ERV$query == genome_index[a,1], ]
	temp_pos <- list()
	for(b in 1:nrow(a_temp)) {
		temp_pos[[b]] <- seq(from=a_temp$q_start[b], to=a_temp$q_end[b], by=1)
	}
	ERV_pos[[a]] <- unlist(temp_pos)


	# get locations of all CDS per scaffold
	a_temp <- cds[cds[,1] == genome_index[a,1], ]
	temp_pos <- list()
	for(b in 1:nrow(a_temp)) {
		temp_pos[[b]] <- seq(from=a_temp[b,4], to=a_temp[b,5], by=1)
	}
	CDS_pos[[a]] <- unlist(temp_pos)
}

# loop for each window to get summarize content
TE_content <- list()
ERV_content <- list()
CR1_content <- list()
CDS_content <- list()
for(a in 1:nrow(window_matrix)) {
	if(a %% 1000 == 0) {print(a)}
	position_index <- match(window_matrix$scaffold[a], genome_index[,1])
	
	# get number of TE sites in this window
	a_temp <- TE_pos[[position_index]]
	a_temp <- a_temp[a_temp >= window_matrix$start[a] & a_temp <= window_matrix$end[a]]
	TE_content[[a]] <- length(a_temp) / window_size
	
	# get number of ERV sites in this window
	a_temp <- ERV_pos[[position_index]]
	a_temp <- a_temp[a_temp >= window_matrix$start[a] & a_temp <= window_matrix$end[a]]
	ERV_content[[a]] <- length(a_temp) / window_size
	
	# get number of CR1 sites in this window
	a_temp <- CR1_pos[[position_index]]
	a_temp <- a_temp[a_temp >= window_matrix$start[a] & a_temp <= window_matrix$end[a]]
	CR1_content[[a]] <- length(a_temp) / window_size
	
	# get number of CDS sites in this window
	a_temp <- CDS_pos[[position_index]]
	a_temp <- a_temp[a_temp >= window_matrix$start[a] & a_temp <= window_matrix$end[a]]
	CDS_content[[a]] <- length(a_temp) / window_size
}

# combine all
output <- data.frame(scaffold=as.character(window_matrix[,1]), start=as.numeric(window_matrix[,2]), end=as.numeric(window_matrix[,3]), cds=as.numeric(unlist(CDS_content)), repeats=as.numeric(unlist(TE_content)), cr1=as.numeric(unlist(CR1_content)), erv=as.numeric(unlist(ERV_content)))

# write output
write.table(output, file="wbnu_window_summary.txt", sep="\t", quote=F, row.names=F)











