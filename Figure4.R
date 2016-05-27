library(ape)
library(adephylo)

DIR <- "Molecular-clock_figures-and-data"

outFile <- paste(DIR, "/results/Figure4.pdf", sep="")
pdf(outFile,width=(6*3), height=(6*1.5))
par(mfrow=c(1,2),oma=c(0,0,3.5,0))
#layout(matrix(c(1,2), nrow = 1), heights=c(0.5))
par(mar=c(5.1,4.1,4.1,2.1))

## human chimpanzee comparison
mytree <- function(col) {
	
cat("((human:",file[1,col],",chimpanzee:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],");", sep = "",file = "tree.tre")
}

# rates
file <- read.table(paste(DIR,"data/Figure4a.txt", sep=""), header = T, sep="\t")
outfile <- "Figure4a"
col <- which(colnames(file)=="Cpg_ts_rate")
mytree(col)
tree <- read.tree("tree.tre")

d_HC = file[1,col]  
G <- 28.4  # Average generation time, Kong et al.
mu <- 1.12e-7/G  # CpG Transition mutation rate, Kong et al.
t_d <- d_HC/(mu)

plot(tree, main=paste("(a) Human-chimpanzee divergence: ",round(t_d/1e6,1), " Mya",sep=""),edge.width = 1.5,cex=2, cex.main=2.4,use.edge.length=TRUE,label.offset=0.002)
labels <- c("",paste(round(tree$edge.length*100,3), "%",sep="")[2:3],"")
edgelabels(labels,adj=c(0.5,-0.5),frame="none", cex=2, offset = 0.2)

# human gorilla comparison
mytree <- function(col) {	
cat("((human:",file[1,col],",gorilla:",file[2,col],"):",file[3,col],",orangutan:",file[4,col],");", sep = "",file = "tree.tre")
}

file <- read.table(paste(DIR,"data/figure4b.txt", sep=""), header = T, sep="\t")
col <- which(colnames(file)=="Cpg_ts_rate")
mytree(col)

d_HC =  file[1,col] 
G <- 28.4
mu <- 1.12e-7/G
t_d <- d_HC/(mu)

tree <- read.tree("tree.tre")
plot(tree, main=paste("(b) Human-gorilla divergence: ",round(t_d/1e6,1), " Mya",sep=""),edge.width = 1.5,cex=2, cex.main=2.4,use.edge.length=TRUE,label.offset=0.002)
labels <- c("",paste(round(tree$edge.length*100,3), "%",sep="")[2:3],"")
edgelabels(labels,adj=c(0.5,-0.5),frame="none", cex=2, offset = 0.2)
dev.off()
