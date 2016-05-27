#Figure S7
library(ape)
library(adephylo)

mytree <- function(col) {
 	 cat("((((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],"):",file[5,col],",(baboon:",file[6,col],",rhesus_macaque:",file[7,col],"):",file[8,col],"):",file[9,col],",marmoset:",file[10,col],");", sep="", file="tree.tre")	
}

DIR <- "Molecular-clock_figures-and-data"
dataset <- c("a","b")

outFile <- paste(DIR, "results/FigureS7.pdf", sep="")
pdf(outFile,width=(6*1.5), height=(6*1.8))
par(mfrow=c(2,1),oma=c(0,0,3.5,0))
layout(matrix(c(1,2), nrow = 2), heights=c(0.5,0.5))

removeSpecies <- c("marmoset") 

for (n in 1:length(dataset)) {
	file <- read.table(paste(DIR,"data/FigureS6-7",dataset[n],".txt", sep=""), header = T, sep="\t")
	results <- NA
	results <- matrix(nrow=6, ncol=2)
	rownames(results) <- c("A/T", "G/C","CpG \noutside CGI", "CpG \nin CGI", "non-CpG G/C\n outside CGI", "non-CpG G/C\n in CGI")
	colnames(results) <- c("transitions","transversions") 

	col <- 1
	for (mrow in 1:nrow(results)) {
		for (mcol in 1:ncol(results)) {
			tree <- NA
			tree_rmSpecies <- NA
			dist_root <- NA
		
			col <- col+1
			mytree(col)
			tree <- read.tree("tree.tre")
			
			tree_rmSpecies <- drop.tip(tree, removeSpecies)
			dist_root <- distRoot(tree_rmSpecies)
			results[mrow,mcol] <- var(dist_root/mean(dist_root))
			#unlink("tree.tre")
		}	
	}
	title <- species[n]
	par(mar = c(7, 4, 2, 2) + 0.9)
	plot(results[,1],col="red",pch=16, ylab="Root-leaf variance",xlab="Substitution type",xaxt="n", ylim=c(0,max(results)))
	points(results[,2],col="blue",pch=16, cex.axis=1,cex.lab=1.2)
#axis(1, at=1:nrow(results), labels=rownames(results), cex.axis=0.7)
	axis(1, at=1:nrow(results), labels=FALSE, tck=-0.02)
	text(seq(1,nrow(results),1),par("usr")[3]-0.005, srt = 0,adj= 0.5, xpd = TRUE,labels=rownames(results),cex=1)
	legend("topright",col=c("red","blue"),pch=16,legend=c("transitions", "transversions"),cex=1)

}
mtext("     (a) Multiz", outer = TRUE, cex =1.5, adj = 0, font=2)
mtext("(b) EPO  ", line= 3.8,cex=1.5, adj = -0.1, font=2)

dev.off()
