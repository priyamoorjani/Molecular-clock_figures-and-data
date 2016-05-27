# Figure 2
library(ape)
library(adephylo)

mytree <- function(col) {
	
cat("((((((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],"):",file[5,col],",(((rhesus_macaque:",file[6,col],",crab-eating_macaque:",file[7,col],"):",file[8,col],",baboon:",file[9,col],"):",file[10,col],",green_monkey:",file[11,col],"):",file[12,col],"):",file[13,col],",(marmoset:",file[14,col],",squirrel_monkey:",file[15,col],"):",file[16,col],"):",file[17,col],",bushbaby:",file[18,col],"):",file[19,col],",mouse:",file[20,col],");", sep="", file="tree.tre")
}

DIR <- "Molecular-clock_figures-and-data"

file <- read.table(paste(DIR,"data/Figure2-3a.txt", sep=""), header = T, sep="\t")
mut_type <- c("Cpg_ts_rate","NoncpgGC_ts_rate")
mut_label <- c("CpG transitions","non-CpG G/C transitions")
species_type <- c("OWM","NWM")

outFile <- paste(DIR, "/results/Figure2.pdf", sep="")
pdf(outFile,width=(6*2.2), height=(6*1.8))
par(mfrow=c(2,2),oma=c(0,0,3.5,0))
layout(matrix(c(1,3,2,4), nrow = 2), heights=c(0.5,0.5))
par(mar=c(5.1,4.1,4.1,2.1))

for (n in 1:length(species_type)) {
	if (species_type[n] == "OWM") { 
		# change branch length color 
		for (i in seq(1,5,1)) { vec.color[i] <- "blueviolet"}
		for (i in seq(6,12,1)) { vec.color[i] <- "lightseagreen"}
		removeSpecies <- c("gibbon","marmoset","squirrel_monkey","bushbaby", "mouse")
	} else {		
		# change branch length color 
		for (i in seq(1,5,1)) { vec.color[i] <- "blueviolet"}
		for (i in seq(6,12,1)) { vec.color[i] <- "coral2"}
		removeSpecies <- c("bushbaby", "mouse", "rhesus_macaque", "crab-eating_macaque","baboon","green_monkey")
	}
	
	for (m in 1:length(mut_type)) {
		col <- which(colnames(file)==mut_type[m])			

		mytree(col)
		tree <- read.tree("tree.tre")
		tree_rmSpecies <- drop.tip(tree, removeSpecies)

		dist_root <- distRoot(tree_rmSpecies)
		ratio <- mean(dist_root[4:length(dist_root)])/ mean(dist_root[1:3])
		title <- paste(mut_label[m])
		norm_edge <- dist_root/dist_root[1]
		node_label <- paste(tree_rmSpecies$tip.label,": ",format(round(norm_edge,2),nsmall=2),sep="")
		node_label[1] <- paste(node_label[1]," ( = ",format(round(dist_root[1],4),nsmall=2),")",sep="")
		tree_rmSpecies$tip.label <-node_label
		plot(tree_rmSpecies, edge.color = vec.color, main=title,edge.width = 1.5,cex=1.6, cex.main=1.7, font.main=1,use.edge.length=TRUE)
		add.scale.bar(cex=1.3)
		unlink("tree.tre")
	}
}
mtext("     (a) Hominoid-OWM", outer = TRUE, cex =1.5, adj = 0, font=2)
mtext("(b) Hominoid-NWM  ", line= 3.8,cex=1.5, adj = -2.2, font=2)
dev.off()

