# Figure S13
DIR <- "Molecular-clock_figures-and-data"

file <- read.table(paste(DIR, "/data/FigureS13.txt",sep=""),header=TRUE, sep="\t")
species <- c("human","chimpanzee","orangutan","rhesus\nmacaque","crab-eating\nmacaque","baboon","green\nmonkey","marmoset","squirrel\nmonkey")
data <- file[,2:ncol(file)]
data <- data/data$Cpg_ts_rate
drop <- c("Cpg_ts_rate")
data <- data[ , !(names(data) %in% drop)]
data <- t(data)
colnames(data) <- species
label <- c("A/T transversions", "A/T transitions",   "non-CpG G/C tranversions", "non-CpG G/C transitions",  "CpG transversions")
colors <- c("red","gold","springgreen3","dodgerblue","hotpink")

outFile <- paste(DIR, "/results/FigureS13.pdf", sep="")
pdf(outFile,width=(6*2.5), height=(6*1.5))
layout(rbind(1,2), heights=c(7,0.5))  # put legend on bottom 1/8th of the chart
par(mar=c(5.1,5.1,4.1,2.1))
barplot(as.matrix(data), xlab="", cex.lab=1.4, col=colors, ylab="Relative mutability of different mutation types", border = FALSE, cex.axis=1.2, cex.names=1.2, beside=TRUE)
par(mar=c(0, 0, 0, 0))
# c(bottom, left, top, right)
plot.new()

legend('center','groups', label, pch=15, col=colors,ncol=3,bty ="n", cex=1.4)
dev.off()