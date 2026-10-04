# Triangle plots of hybrid index vs interclass heterozygosity for P. alecto and
# P. conspicillatus, based on ancestry-informative markers (AIMs) with an allele
# frequency difference >= 0.85 between the species (triangulaR).
#
# Run from this directory:  cd 3_SNP_Structure/scripts && Rscript triangle_hybridization_plot.R
#
# Inputs:
#   ../../1_SNP_Calling/data/palecto_pconspicillatus_narrow.vcf.gz  "alecto ex alecto alecto + conspicillatus" SNP set
#   ../../0_Metadata/data/palecto_pconspicillatus.popmap
# Outputs:
#   ../figures/Triangle_plot_BFF_SFF.pdf              Figure S10
#   ../figures/Triangle_plot_BFF_SFF_missingness.pdf  Figure S11
# Outlier labels were added in a graphics editor.

library(ggplot2)
library(dplyr)
library(triangulaR)
library(vcfR)

bffsff_vcf <- "../../1_SNP_Calling/data/palecto_pconspicillatus_narrow.vcf.gz"
data <- read.vcfR(bffsff_vcf, verbose = F)
bffsff_popfile <- read.csv("../../0_Metadata/data/palecto_pconspicillatus.popmap", header = FALSE, sep = "\t")
colnames(bffsff_popfile) <- c("id","pop")
# modify popfile to remove outgroups and merge species
bffsff_popfile_sp <- bffsff_popfile %>% filter(pop!="INDO_outgroup") %>% mutate(pop=ifelse(pop=="Wet_tropics"| pop=="PNG","SFF","BFF"))

# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
bffsff.diff85 <- alleleFreqDiff(vcfR = data, pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF", difference = 0.85)
hi.het85 <- hybridIndex(vcfR = bffsff.diff85,  pm = bffsff_popfile_sp, p1 = "BFF", p2 = "SFF")
# Exclude individuals missing one third or more of the AIMs (>= 5 of 15)
hi.het85 <- hi.het85 %>% filter(perc.missing < 0.33)

dir.create("../figures", showWarnings = FALSE)
triangle.plot(hi.het85,colors = c("#E495A5","#C29DDE"),ind.labels = T,alpha=0.75)
ggsave("../figures/Triangle_plot_BFF_SFF.pdf",width = 5,height=4)
missing.plot(hi.het85)
ggsave("../figures/Triangle_plot_BFF_SFF_missingness.pdf",width = 5.5,height=4)
