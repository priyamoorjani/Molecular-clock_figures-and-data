# Figure S1

DIR <- "Molecular-clock_figures-and-data"

outFileName <- paste("results/FigureS1.png", sep='')

for (chr in seq(1,22,1)) {
        fileName <- paste(DIR,"data/GSE30340_classify_cpgIsland_chr", chr, ".txt.gz", sep="")
        file <- read.table(fileName, header = F)
        if (chr == 1) {
                inIsland <- file[ which(file$V6=="In_CpG_Island"), ]
                notIsland <- file[ which(file$V6=="Not_in_CpG_Island"), ]
        } else {
                inIsland <- rbind(inIsland,file[ which(file$V6=="In_CpG_Island"), ])
                notIsland <- rbind(notIsland,file[ which(file$V6=="Not_in_CpG_Island"), ])
        }
}

png(outFileName, width=900, height=450)
par(mfrow=c(1,2))
title <- paste('CpG sites in CGI',sep='')
hist(inIsland$V4, col="salmon",main=title, xlab="methylation level", ylim=c(0,4e6))
title <- paste('CpG sites outside CGI',sep='')
hist(notIsland$V4, col="dodger blue", main=title, xlab="methylation level")
dev.off()