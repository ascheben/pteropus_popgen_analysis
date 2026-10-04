# Plot fastSTRUCTURE admixture proportions (K = 4) for P. alecto and
# P. conspicillatus, grouped by geographic region.
#
# Run from this directory:  cd 3_SNP_Structure/scripts && Rscript structure_plot.R
#
# Inputs:
#   ../data/structure/structure_palecto_pconspicillatus_out.<K>.meanQ   fastSTRUCTURE output, K = 1-15
#   ../data/structure/structure_palecto_pconspicillatus_labels.txt      sample labels (row order of .meanQ)
#   ../data/structure/*.fam                    PLINK sample IDs (row order of .meanQ) used to match the popmap
#   ../../0_Metadata/data/palecto_pconspicillatus.popmap
# Output:
#   ../figures/Admixture_barplot_faststructure_bff_sff_k4.pdf   Figure 3E
#   (region and species labels were finalised in a graphics editor)

library(pophelper)
library(gridExtra)

inpath <- "../data/structure/"
labels <- "../data/structure/structure_palecto_pconspicillatus_labels.txt"
famfile <- list.files(path = inpath, pattern = "\\.fam$", full.names = TRUE)
grplabels <- "../../0_Metadata/data/palecto_pconspicillatus.popmap"

shiny <- c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E")

#Add all files in directory to list
sfiles <- list.files(path = inpath, pattern = "structure_palecto_pconspicillatus_out.*.meanQ$", full.names=TRUE)

#read files in from list
slist <- readQ(files=sfiles,indlabfromfile=T)
#read individual labels
inds <- read.delim(labels,header=F,stringsAsFactors=FALSE)
#add ind names as rownames to all tables
if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"rownames<-",inds$V1)

# Assign regions using the PLINK sample IDs (same order as the labels)
inds$sample <- read.table(famfile, stringsAsFactors = FALSE)$V2
grouplabset <- read.delim(grplabels, header=F,stringsAsFactors=F)
colnames(grouplabset) <-  c("name","region")
inds$id  <- 1:nrow(inds)
shared_labs <- merge(inds,grouplabset,by.x = c("sample"), by.y=c("name"),all.x = TRUE)
shared_labs <- shared_labs[order(shared_labs$id), ]

regions <- as.data.frame(shared_labs[, c("region")])
# K = 4 was selected with fastSTRUCTURE chooseK.py
k4 <- grep("_out\\.4\\.meanQ$", sfiles)
p1 <- plotQ(slist[k4],returnplot=T,exportplot=F,basesize=11,
            linesize=0.8,
            pointsize=3,
            showindlab=T,
            useindlab = T,
            grplab =regions,
            sortind = "all",
            indlabsize=4,
            ordergrp=T,
            grplabsize=3.5,
            clustercol=shiny)

dir.create("../figures", showWarnings = FALSE)
pdf("../figures/Admixture_barplot_faststructure_bff_sff_k4.pdf", width = 12, height = 5);
grid.arrange(p1$plot[[1]])
dev.off()
