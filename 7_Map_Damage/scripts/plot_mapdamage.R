# Compare DNA damage (C to T substitution at the first 5' read position after
# clipping) between extant (CSIRO, n = 66) and historical (WA museum, n = 15)
# P. alecto samples.
#
# Run from this directory:  cd 7_Map_Damage/scripts && Rscript plot_mapdamage.R
#
# Inputs:
#   ../data/5pCtoT_freq_clipped_all.txt          mapDamage 5pCtoT_freq.txt of all samples (pos, CtoT_fraction, sample)
#   ../../0_Metadata/data/pteropus_metadata.txt  sample metadata (Species, Source)
# Output:
#   ../figures/CtoT_extant_vs_museum.pdf         Figure S2 (group labels renamed to Extant/Historical in the manuscript)

library(ggplot2)
library(dplyr)
library(ggpubr)

pteropus_meta <-read.csv("../../0_Metadata/data/pteropus_metadata.txt",  header = TRUE, sep = "\t")
ctot <- read.csv("../data/5pCtoT_freq_clipped_all.txt.gz",header =TRUE, sep = "\t")

# Reads were clipped by 5 bp, so position 6 of the mapDamage output is the first
# retained 5' base (position 1 of the clipped reads)
ctot_pos6 <- ctot[ctot$pos==6,]
ctot_pos6 <- merge(ctot_pos6,pteropus_meta,by.x = "sample", by.y="Sample.Identifier")
ctot_pos6_bff <- ctot_pos6[ctot_pos6$Species=="BFF",]

# Figure S2 compares CSIRO and WA museum samples; the 5 bat carer samples are not shown
ctot_pos6_bff %>% filter(Source %in% c("CSIRO","WA museum")) %>% ggplot(aes(Source,CtoT_fraction)) +
  ylab("C to T conversion fraction (5' position 1)")+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_bw()+
  stat_compare_means(method="t.test",label.y=0.0075)

dir.create("../figures", showWarnings = FALSE)
ggsave("../figures/CtoT_extant_vs_museum.pdf")
