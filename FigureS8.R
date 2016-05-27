# Figure S8
library(ape)
library(adephylo)

mytree <- function(col) {
 	 cat("((((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],"):",file[5,col],",(baboon:",file[6,col],",rhesus_macaque:",file[7,col],"):",file[8,col],"):",file[9,col],",marmoset:",file[10,col],");", sep="", file="tree.tre")	
}

DIR <- "Molecular-clock_figures-and-data"
dataset <- c("a","b")

removeSpecies <- c("marmoset")
title <- c("W->W","S->S",	"CpG\n S->S", "non-CpG G/C\n S-> S",	"W->S","S->W","CpG\n S->W","non-CpG G/C\n S->W")
popSymbol <- c(rep(15,4), rep(17,4))

outFile <- paste(DIR, "results/FigureS8.pdf", sep="")
pdf(outFile,width=(6*1.5), height=(6*1.8))
par(mfrow=c(2,1),oma=c(0,0,3.5,0))
layout(matrix(c(1,2), nrow = 2), heights=c(0.5,0.5))

for (n in 1:length(dataset)) {
	results <- NA; file1 <- NA; file2 <- NA
	results <- matrix(nrow=length(title),ncol=1)
	rownames(results) <- title
	colnames(results) <- c("RLV")

	file <- read.table(paste(DIR,"data/FigureS8",dataset[n],".txt", sep=""), header = T, sep="\t")

	col <- 1
	for (mcol in 1:ncol(results)) {	
		for (mrow in 1:nrow(results)) {
			col <- col + 1
			tree <- NA
			tree_rmSpecies <- NA
			dist_root <- NA
		
			mytree(col)
			tree <- read.tree("tree.tre")
			tree_rmSpecies <- drop.tip(tree, removeSpecies)
			dist_root <- distRoot(tree_rmSpecies)
			results[mrow,mcol] <- var(dist_root/mean(dist_root))
			unlink("tree.tre")	
		}	
	}
	par(mar = c(7, 4, 2, 2)+0.9)
	plot(1:nrow(results), results[,1],col=c("burlywood4"),pch=popSymbol, ylab="Root-leaf variance",xlab="Substitution type",xaxt="n",ylim=c(0,max(results)+0.01),cex.lab=1.1,cex.axis=1)
	axis(1, at=1:nrow(results), labels=FALSE, tck=-0.01)
	text(seq(1,nrow(results),1),par("usr")[3]-0.004, srt = 0,adj= 0.5, xpd = TRUE,labels=rownames(results),cex=1)
	legend("topright",col=c("burlywood4","burlywood4"),pch=c(15,17),legend=c("Not sensitive to BGC","Sensitive to BGC" ),cex=1)
}
mtext("     (a) Multiz", outer = TRUE, cex =1.5, adj = 0, font=2)
mtext("(b) EPO  ", line= 3.8,cex=1.5, adj = -0.1, font=2)

dev.off()

