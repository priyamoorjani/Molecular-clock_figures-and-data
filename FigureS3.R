# Figure S3
library(ape)
library(adephylo)
DIR <- "Molecular-clock_figures-and-data"
data <- c("Figure2-3a","FigureS3-4_G","FigureS3-4_GI")
list <- c("noGnoGI", "noCnoGI", "noGnoO")

mytree <- function(col, indicator) {		
switch(indicator,
noGnoGI=  { cat("((((((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],"):",file[5,col],",(((rhesus_macaque:",file[6,col],",crab-eating_macaque:",file[7,col],"):",file[8,col],",baboon:",file[9,col],"):",file[10,col],",green_monkey:",file[11,col],"):",file[12,col],"):",file[13,col],",(marmoset:",file[14,col],",squirrel_monkey:",file[15,col],"):",file[16,col],"):",file[17,col],",bushbaby:",file[18,col],"):",file[19,col],",mouse:",file[20,col],");", sep="", file="tree.tre") },

noCnoGI=  { cat("((((((human:",file[1,col],",gorilla:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],"):",file[5,col],",(((rhesus_macaque:",file[6,col],",crab-eating_macaque:",file[7,col],"):",file[8,col],",baboon:",file[9,col],"):",file[10,col],",green_monkey:",file[11,col],"):",file[12,col],"):",file[13,col],",(marmoset:",file[14,col],",squirrel_monkey:",file[15,col],"):",file[16,col],"):",file[17,col],",bushbaby:",file[18,col],"):",file[19,col],",mouse:",file[20,col],");", sep="", file="tree.tre") },

noGnoO=
{ cat("((((((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",gibbon:",file[4,col],"):",file[5,col],",(((rhesus_macaque:",file[6,col],",crab-eating_macaque:",file[7,col],"):",file[8,col],",baboon:",file[9,col],"):",file[10,col],",green_monkey:",file[11,col],"):",file[12,col],"):",file[13,col],",(marmoset:",file[14,col],",squirrel_monkey:",file[15,col],"):",file[16,col],"):",file[17,col],",bushbaby:",file[18,col],"):",file[19,col],",mouse:",file[20,col],");", sep="", file="tree.tre") }
)
} 

outFile <- paste(DIR, "/results/FigureS3.pdf", sep="")
pdf(outFile,width=(6*3), height=(6*1.8))
par(mfrow=c(2,2),oma=c(0,0,3.5,0))
layout(matrix(c(1,2,3,4,5,6), nrow = 2), heights=c(0.5,0.5))
par(mar=c(5.1,4.1,4.1,2.1))

removeSpecies <- c("marmoset","squirrel_monkey","bushbaby", "mouse")

for (i in 1:length(list)) {
	species <- list[i]
	file <- read.table(paste(DIR,"data/",data[i],".txt", sep=""), header = T, sep="\t")
	
	mut_type <- c("Cpg_ts_rate","NoncpgGC_ts_rate")
	# cpg sites
	for (m in 1:length(mut_type)) {
		col <- which(colnames(file)==mut_type[m])
		mytree(col, species)
		tree <- read.tree("tree.tre")
		tree_rmSpecies <- drop.tip(tree, removeSpecies)

		# change color of desired edges. Here color each edge to bushbaby in red.
		vec.color <- rep(1,length=nrow(tree_rmSpecies$edge))
		for (i in seq(1,5,1)) { vec.color[i] <- "blueviolet"}
		for (i in seq(6,12,1)) { vec.color[i] <- "lightseagreen"}

		dist_root <- distRoot(tree_rmSpecies)	
		ratio <- mean(dist_root[4:length(dist_root)])/ mean(dist_root[1:3])
		norm_edge <- dist_root/dist_root[1]
		node_label <- paste(tree_rmSpecies$tip.label,": ",format(round(norm_edge,2),nsmall=2),sep="")
		node_label[1] <- paste(node_label[1]," ( = ",format(round(dist_root[1],4),nsmall=2),")",sep="")
		tree_rmSpecies$tip.label <-node_label
		title <- paste("root_OWM / root_hominoid: ",format(round(ratio,2), nsamll=2),sep="")
		plot(tree_rmSpecies, edge.color = vec.color, main=title,edge.width = 1.5,cex=2, cex.main=2, 	font.main=1,use.edge.length=TRUE)
		add.scale.bar(cex=2)
	
		}
}	
mtext("     (a) Transitions at CpG sites", outer = TRUE, cex =1.5, adj = 0, font=2)
mtext("(b) Transitions at non-CpG G/C sites", line= 5,cex=1.5, adj = -15.0, font=2)
dev.off()