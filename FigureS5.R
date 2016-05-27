# Figure S5
 library(ape)
 library(adephylo)
 
mytree <- function(col) {
	
 cat("((((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],"):",file[5,col],",(baboon:",file[6,col],",rhesus_macaque:",file[7,col],"):",file[8,col],"):",file[9,col],",marmoset:",file[10,col],");", sep="", file="tree.tre")
 }
 
DIR <- "Molecular-clock_figures-and-data"
file <- read.table(paste(DIR,"data/FigureS5.txt", sep=""), header = T, sep="\t")
for (i in 1:ncol(file)) {
	if (colnames(file)[i] %in% "average_rate") {
			col <- i
	}
}
mytree(col)
tree <- read.tree("tree.tre")

# set each edge to black color = 1
vec.color <- rep(1,length=nrow(tree$edge))

pdf(paste(DIR, "/results/FigureS5.pdf", sep=""))
plot(tree, edge.color = vec.color, tip.color=c("blueviolet","blueviolet","blueviolet","lightseagreen","lightseagreen","coral2","coral2"), use.edge.length=TRUE,label.offset=0.0005)
legend("bottomleft",col=c("blueviolet","lightseagreen","coral2","brown3"),lty=1,legend=c("Apes","Old World Monkeys","New World Monkeys"),cex=1)
add.scale.bar(0.001,2,cex=1)
dev.off()