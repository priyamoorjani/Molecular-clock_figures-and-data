# Figure S2
# Author: Carlos EG Amorim
DIR <- "Molecular-clock_figures-and-data"

rec <- read.table(paste(DIR, "/data/FigureS2.txt",sep=""),header=TRUE)
both<-rbind(rec$CpG,rec$GC)

pdf(paste(DIR, "/results/FigureS2.pdf", sep=""))
barplot(both,beside=TRUE,xlab="Recombination rate", main="Proportion of G/C sites",col=c("darkgoldenrod1","turquoise3"), legend = c("CpG","nonCpG G/C"), names.arg=c("0", "(0-0.01]", "(0.01-0.1]", "(0.1-1]", "(1-10]", "(10-100]", ">100"))
dev.off()