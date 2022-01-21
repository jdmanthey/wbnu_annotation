options(scipen=999)

# read in window summary
x <- read.table("wbnu_window_summary.txt", header=T, stringsAsFactors=F)

# define window size
window_size <- x$end[1] - x$start[1] + 1

# define number of windows for line plots
num_windows <- 10

# set up row numbers for line plots
line_rows <- seq(from=1,to=nrow(x), by=1)[seq(from=1,to=nrow(x), by=1) %% num_windows == 0]

# what are the unique chromosomes and their bounding areas for plotting?
total_windows <- nrow(x)
scaffold <- unique(x$scaffold)
chr_polygons <- list()
# make the plotting polygons
for(a in 1:length(scaffold)) {
	a1 <- as.numeric(rownames(x))[x[,1] == scaffold[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 0.5), c(a1, 0.5), c(a1, 0))
}

########################################################################
########################################################################
########################################################################
########################################################################
# set up plotting dimensions
par(mfrow=c(2,1))
par(mar=c(1,5,1,0))

########################################################################
########################################################################
########################################################################
########################################################################
# plot CDS
plot(c(-1,-1), ylim=c(0,0.2), xlim=c(1, total_windows), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="CDS Content")
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
# plot
points(as.numeric(rownames(x)), x$cds, pch=19, cex=0.1, col="gray85")

# plot sliding mean line plots
x_windows <- as.numeric(rownames(x))
line_x_axis <- line_rows
line_y_axis <- list()
line_scaffold <- list()
for(b in line_rows) {
	line_y_axis[[b]] <- mean(na.omit(x$cds[x_windows %in% (b - (floor(num_windows / 2))):(b + (floor(num_windows / 2))) & x$scaffold == x$scaffold[as.numeric(rownames(x)) == b]][1:num_windows]))
	line_scaffold[[b]] <- x$scaffold[x_windows == b] 
}
line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
# plot each scaffold at a time (so lines don't connect between scaffolds)
for(a in 1:length(unique(scaffold))) {
	a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[a],]
	lines(a_rep[,1:2], lwd=0.8, col="gray55")
}
########################################################################
########################################################################
########################################################################
########################################################################
# plot repeats
plot(c(-1,-1), ylim=c(0,0.4), xlim=c(1, total_windows), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Repeat Content")
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
# plot
points(as.numeric(rownames(x)), x$repeats, pch=19, cex=0.1, col="gray85")

# plot sliding mean line plots
x_windows <- as.numeric(rownames(x))
line_x_axis <- line_rows
line_y_axis <- list()
line_scaffold <- list()
for(b in line_rows) {
	line_y_axis[[b]] <- mean(na.omit(x$repeats[x_windows %in% (b - (floor(num_windows / 2))):(b + (floor(num_windows / 2))) & x$scaffold == x$scaffold[as.numeric(rownames(x)) == b]][1:num_windows]))
	line_scaffold[[b]] <- x$scaffold[x_windows == b] 
}
line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
# plot each scaffold at a time (so lines don't connect between scaffolds)
for(a in 1:length(unique(scaffold))) {
	a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[a],]
	lines(a_rep[,1:2], lwd=0.8, col="gray55")
}







