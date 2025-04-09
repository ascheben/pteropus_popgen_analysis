library(ggrepel)
library(reshape2)
library(ggplot2)
library(dplyr)
library(colorspace)
library(cowplot)
library(triangulaR)
library(vcfR)


bffsff_vcf <- "../../1_SNP_Calling/data/palecto_pconspicillatus_narrow.vcf.gz"
data <- read.vcfR(bffsff_vcf, verbose = F)
bffsff_popfile <- read.csv("../../0_Metadata/data/palecto_pconspicillatus.popmap", header = FALSE, sep = "\t")
colnames(bffsff_popfile) <- c("id","pop")
# modify popfile to remove outgroups and merge species
bffsff_popfile_sp <- bffsff_popfile %>% filter(pop!="INDO_outgroup") %>% mutate(pop=ifelse(pop=="Wet_tropics"| pop=="PNG","SFF","BFF"))

# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
#bffsff.diff75 <- alleleFreqDiff(vcfR = data, pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF", difference = 0.75)
#bffsff.diff8 <- alleleFreqDiff(vcfR = data, pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF", difference = 0.8)
bffsff.diff85 <- alleleFreqDiff(vcfR = data, pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF", difference = 0.85)
#bffsff.diff9 <- alleleFreqDiff(vcfR = data, pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF", difference = 0.9)
#bffsff.diff95 <- alleleFreqDiff(vcfR = data, pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF", difference = 0.95)


#hi.het75 <- hybridIndex(vcfR = bffsff.diff75,  pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF")
#hi.het8 <- hybridIndex(vcfR = bffsff.diff8,  pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF")
hi.het85 <- hybridIndex(vcfR = bffsff.diff85,  pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF")
#hi.het9 <- hybridIndex(vcfR = bffsff.diff9,  pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF")
#hi.het95 <- hybridIndex(vcfR = bffsff.diff95,  pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF")

#triangle.plot(hi.het75)
#triangle.plot(hi.het8)
#triangle.plot(hi.het85)
#triangle.plot(hi.het9)
#triangle.plot(hi.het95)

# final
triangle.plot(hi.het85,colors = c("#E495A5","#C29DDE"),ind.labels = T,alpha=0.75)
ggsave("../figures/Triangle_plot_BFF_SFF.pdf",width = 5,height=4)
missing.plot(hi.het85)
ggsave("../figures/Triangle_plot_BFF_SFF_missingness.pdf",width = 5.5,height=4)


