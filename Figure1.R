# Figure 1 - overall substitution rates in Multiz
library(ape)
library(adephylo)

mytree <- function(col) {
	
cat("((((((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],"):",file[5,col],",(((rhesus_macaque:",file[6,col],",crab-eating_macaque:",file[7,col],"):",file[8,col],",baboon:",file[9,col],"):",file[10,col],",green_monkey:",file[11,col],"):",file[12,col],"):",file[13,col],",(marmoset:",file[14,col],",squirrel_monkey:",file[15,col],"):",file[16,col],"):",file[17,col],",bushbaby:",file[18,col],"):",file[19,col],",mouse:",file[20,col],");", sep="", file="tree.tre")
}

DIR <- "Molecular-clock_figures-and-data"

file <- read.table(paste(DIR,"data/Figure1.txt", sep=""), header = T, sep="\t")

col <- which(colnames(file)=="average_rate")
mytree(col)
tree <- read.tree("tree.tre")

row.names(tree$edge) <- seq(1,nrow(tree$edge),1)

# set each edge to black color = 1
vec.color <- rep(1,length=nrow(tree$edge))

# Set outgroup color
vec.color[1] <- "gray88"
vec.color[nrow(tree$edge)] <- "gray88"

pdf(paste(DIR, "/results/Figure1.pdf", sep=""))
plot(tree, edge.color = vec.color, tip.color=c("blueviolet","blueviolet","blueviolet","lightseagreen","lightseagreen","lightseagreen","lightseagreen","coral2","coral2","brown3","black"), use.edge.length=TRUE,label.offset=0.002)
legend("bottomleft",col=c("blueviolet","lightseagreen","coral2","brown3","black"),lty=1,legend=c("Apes","Old World Monkeys","New World Monkeys","Prosimians","Outgroup"),cex=1)
add.scale.bar(0.001,3.5,cex=1)
dev.off()
