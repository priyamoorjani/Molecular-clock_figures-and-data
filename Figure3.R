# Root-tip variance - rates
library(ape)
library(adephylo)

mytree <- function(col) {
	
cat("((((((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],"):",file[5,col],",(((rhesus_macaque:",file[6,col],",crab-eating_macaque:",file[7,col],"):",file[8,col],",baboon:",file[9,col],"):",file[10,col],",green_monkey:",file[11,col],"):",file[12,col],"):",file[13,col],",(marmoset:",file[14,col],",squirrel_monkey:",file[15,col],"):",file[16,col],"):",file[17,col],",bushbaby:",file[18,col],"):",file[19,col],",mouse:",file[20,col],");", sep="", file="tree.tre")
}

DIR <- "Molecular-clock_figures-and-data"
outFile <- paste(DIR, "results/Figure3.pdf", sep="")
pdf(outFile,width=(6*1.5), height=(6*1.8))
par(mfrow=c(2,1),oma=c(0,0,3.5,0))
layout(matrix(c(1,2), nrow = 2), heights=c(0.5,0.5))

# root-leaf variance by mutation type
file <- read.table(paste(DIR,"data/figure2-3a.txt", sep=""), header = T, sep="\t")
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

		removeSpecies <- c("mouse","bushbaby") 
		tree_rmSpecies <- drop.tip(tree, removeSpecies)
		dist_root <- distRoot(tree_rmSpecies)
		results[mrow,mcol] <- var(dist_root/mean(dist_root))
		unlink("tree.tre")
	}	
}

par(mar = c(7, 4, 2, 2) + 0.9)
plot(results[,1],col="red",pch=16, ylab="Root-leaf variance",xlab="Substitution type",xaxt="n", ylim=c(0,max(results)), cex.lab=1.2,cex.axis=1.1)
points(results[,2],col="blue",pch=16)
axis(1, at=1:nrow(results), labels=FALSE, tck=-0.02)
text(seq(1,nrow(results),1),par("usr")[3]-0.007, srt = 0,adj= 0.5, xpd = TRUE,labels=rownames(results),cex=1)
legend("topright",col=c("red","blue"),pch=16,legend=c("transitions", "transversions"),cex=1)


# root-leaf variance- for types subject and not subject to BGC
file <- read.table(paste(DIR,"data/figure3b.txt", sep=""), header = T, sep="\t")
title <- c("W->W","S->S",	"CpG\n S->S", "non-CpG G/C\n S-> S",	"W->S","S->W","CpG\n S->W","non-CpG G/C\n S->W")

results <- NA
results <- matrix(nrow=length(title),ncol=1)
rownames(results) <- title
colnames(results) <- c("RLV")
popSymbol <- c(rep(15,4), rep(17,4))

col <- 1
for (mcol in 1:ncol(results)) {	
	for (mrow in 1:nrow(results)) {
		col <- col +1 
		
		tree <- NA
		tree_rmSpecies <- NA
		dist_root <- NA
		
		mytree(col)
		tree <- read.tree("tree.tre")

		removeSpecies <- c("mouse","bushbaby")

		tree_rmSpecies <- drop.tip(tree, removeSpecies)
		dist_root <- distRoot(tree_rmSpecies)
		results[mrow,mcol] <- var(dist_root/mean(dist_root))
		unlink("tree.tre")	
	}	
}

par(mar = c(7, 4, 2, 2)+0.9)
plot(1:nrow(results), results[,1],col=c("burlywood4"),pch=popSymbol, ylab="Root-leaf variance",xlab="Substitution type",xaxt="n",ylim=c(0,max(results)), cex.lab=1.2,cex.axis=1.1)
axis(1, at=1:nrow(results), labels=FALSE, tck=-0.02)
text(seq(1,nrow(results),1),par("usr")[3]-0.007, srt = 0,adj= 0.5, xpd = TRUE,labels=rownames(results),cex=1)
#axis(1, at=1:nrow(results), labels=rownames(results), cex.axis=0.9, tck=-0.02)
legend("topright",col=c("burlywood4","burlywood4"),pch=c(15,17),legend=c("Not sensitive to BGC","Sensitive to BGC" ),cex=1)

mtext("     (a) variation in substitution rates, by mutation type and context", outer = TRUE, cex =1.3, adj = 0, font=2)
mtext("(b) variation in substitution rates, for types subject and not subject to BGC", line= 3.8,cex=1.3, adj = 1.5, font=2)
dev.off()
