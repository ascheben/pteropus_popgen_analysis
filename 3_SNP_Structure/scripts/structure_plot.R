library(pophelper)
library(gridExtra)
library(dplyr)

inpath <- "../data/structure/"
labels <- "../data/structure/bffsff_labels.txt"
grplabels <- "../../0_Metadata/data/palecto_pconspicillatus.popmap"


lenK <- 15

shiny <- c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E")


#Add all files in directory to list
sfiles <- list.files(path = inpath, pattern = "structure_palecto_pconspicillatus_out.*.meanQ$", full.names=TRUE)

#read files in from list and add individual labes from file
#files can be modifief to show cultivar names
slist <- readQ(files=sfiles,indlabfromfile=T)
#read individual labels
inds <- read.delim(labels,header=F,stringsAsFactors=FALSE)
#add ind names as rownames to all tables
if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"rownames<-",inds$V1)

grouplabset <- read.delim(grplabels, header=F,stringsAsFactors=F)
colnames(grouplabset) <-  c("name","region")
#grouplabset$species <- gsub("SFF","cSFF",grouplabset$species)
inds$id  <- 1:nrow(inds)
shared_labs <- merge(inds,grouplabset,by.x = c("V1"), by.y=c("name"),all.x = TRUE)
shared_labs <- shared_labs[order(shared_labs$id), ]


#p1 <- plotQ(slist[10],returnplot=T,exportplot=F,basesize=9,
#            linesize=0.8,pointsize=3,showindlab=T,useindlab = T,sortind = "all",clustercol=shiny)
#grid.arrange(p1$plot[[1]])

regions <- as.data.frame(shared_labs[, c("region")])
p1 <- plotQ(slist[10],returnplot=T,exportplot=F,basesize=11,
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

pdf("../figures/Admixture_barplot_faststructure_bff_sff_k4.pdf", width = 12, height = 5);    
grid.arrange(p1$plot[[1]])
dev.off()
