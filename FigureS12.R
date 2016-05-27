# Figure S12 Maximum Likelihood - HC
library(ape)
library(adephylo)

mytree <- function(col) {
	
cat("((human:",file[1,col],",gorilla:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],");", sep = "",file = "tree.tre")
}

DIR <- "Molecular-clock_figures-and-data"
file <- read.table(paste(DIR,"data/FigureS12.txt", sep=""), header = TRUE,sep="\t")
title <- c("ancestral A/T","ancestral A/T", "ancestral CpG","ancestral CpG","ancestral non-CpG G/C","ancestral non-CpG G/C")

pdf(paste(DIR, "/results/FigureS12.pdf", sep=""),width=(6*3), height=(6*2))
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(2,2),oma=c(0,0,2,2))
layout(matrix(c(1,2,3,4,5,6), nrow = 2), heights=c(0.4,0.4))
par(mar=c(5.1,4.1,4.1,2.1))

for (col in 2:(ncol(file)))
{
	type <- colnames(file)[col]
	mytree(col)
	tree <- read.tree("tree.tre")
	index <- col-1
	removeSpecies <- c("orangutan")
	tree_rmSpecies <- drop.tip(tree, removeSpecies)
	dist_root <- distRoot(tree_rmSpecies)
	ratio <- dist_root[2]/dist_root[1]
	string <- paste(title[index],": ",format(round(ratio,2),nsmall=2),sep="")
	plot (tree, main = string, label.offset=0.0002, edge.width = 2, cex =2, use.edge.length=TRUE, cex.main=2, font.main=1)
	labels <- paste(format(round(tree$edge.length*100,3),nsmall=3), "%",sep="")
	edgelabels(labels,adj=c(0.5,-0.2),frame="none", cex=1.8, offset = 0.2)	
	unlink("tree.tre")
}
mtext("   (a) Transitions", outer = TRUE, cex =1.5, adj = 0, font=2)
mtext("(b) Transversions ", line= 4,cex=1.5, adj = -4.1, font=2)
dev.off()

