# Plot and test genome-wide diversity (heterozygosity, nucleotide diversity) and
# genetic differentiation (FST) for the four flying-fox species and the regional
# populations of P. alecto and P. conspicillatus.
#
# Run from this directory:  cd 2_SNP_Stats/scripts && Rscript plot_stats.R
#
# Inputs (../data/ unless stated):
#   pteropus_heterozygosity.txt            per-sample heterozygosity (Pop, Sample, Het_frac, ...)
#   pteropus_nucleotide_diversity.txt.gz   pixy pi in 100 kb windows per species and region
#   pteropus.fst_summary.tsv               stacks populations FST between species
#   palecto_pconspicillatus.fst_summary.tsv  stacks populations FST between regions
#   ../../0_Metadata/data/palecto_pconspicillatus.popmap
# Outputs (../figures/):
#   Heterozygosity_species_populations.pdf   Figure S3
#   Nucleotide_diversity_all_species.pdf     Figure S4
#   Fst_species.pdf                          Figure S5
#   Fst_bff_sff_regional.pdf                 Figure S6
#   *_wilcoxon.tsv                      pairwise Wilcoxon rank-sum tests (Tables S4, S6)
#   *_summary.tsv                            mean/sd per group (Tables S3, S5)
# Taxon and region labels were renamed to full names and panels arranged in a
# graphics editor for the manuscript.

library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)

get_box_stats <- function(y, upper_limit = 0.02) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 4), "\n",
      "Median =", round(median(y), 4), "\n"
    )
  ))
}

# Pairwise two-sided Wilcoxon rank-sum tests (Bonferroni-adjusted), as a long table
pairwise_wilcox_table <- function(x, g) {
  p <- pairwise.wilcox.test(x = x, g = g, p.adjust.method = "bonferroni", exact = FALSE)$p.value
  p <- melt(p, na.rm = TRUE)
  colnames(p) <- c("group1", "group2", "p_adj")
  p
}

dir.create("../figures", showWarnings = FALSE)

hetstat <- read.delim("../data/pteropus_heterozygosity.txt",sep = '\t',header=TRUE)
bffsff_pop <- read.delim("../../0_Metadata/data/palecto_pconspicillatus.popmap",sep = '\t',header=FALSE)
# Read Fst tables
fst_species <-  read.delim("../data/pteropus.fst_summary.tsv")
fst_bffsff <-  read.delim("../data/palecto_pconspicillatus.fst_summary.tsv")

#########################
## Nucleotide diversity #
#########################

pi <- read.delim("../data/pteropus_nucleotide_diversity.txt.gz",header = TRUE, sep = "\t")
pi <- pi %>% mutate(group=ifelse(pop == "BFF" | pop == "GHFF" | pop == "LRFF" | pop == "alecto_alecto" | pop == "SFF","species","population")) %>%
  mutate(species=ifelse(pop=="PNG" | pop=="Wet_tropics","SFF","BFF"))
# Remove windows with fewer than 50 sites
pi <- pi %>% filter(no_sites>=50)

pi$pop <- factor(pi$pop, levels = c("alecto_alecto","BFF","SFF","LRFF","GHFF", "PNG","Wet_tropics","INDO","N_AUS","NQ","E_COAST"))
p1 <- ggplot(pi[pi$group=="population",],aes(pop,avg_pi,fill=species)) + geom_boxplot(outliers = FALSE,width=0.75) + theme_bw()+
  scale_fill_manual(values = c("lightblue","firebrick"))+
  xlab("")+
  ylab("π")+
  theme(text = element_text(size=18,angle = 45),axis.text.x = element_text(angle = 45, vjust = 0.5))+
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 2)+
  stat_summary(fun = "mean", geom = "point", shape = 2, size = 2, color = "black")

p2 <- ggplot(pi[pi$group=="species",],aes(reorder(pop,avg_pi),avg_pi,fill="grey")) +
  geom_boxplot(outliers = FALSE,width=0.61) + theme_bw()+
  scale_fill_manual(values = c("grey"))+
  xlab("")+
  ylab("π")+
  theme(text = element_text(size=18),axis.text.x = element_text(angle = 45, vjust = 0.5))+
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 2)+
  stat_summary(fun = "mean", geom = "point", shape = 2, size = 2, color = "black")

plot_grid(p1,p2)
ggsave("../figures/Nucleotide_diversity_all_species.pdf",height = 6,width = 11)

pi_summary <- pi %>% group_by(group, pop) %>% summarize(
  n_windows = n(),
  mean = mean(avg_pi, na.rm = TRUE),
  sd = sd(avg_pi, na.rm = TRUE), .groups = "drop")
write.table(pi_summary, "../figures/Nucleotide_diversity_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Pairwise tests between species and between regional populations
pi_tests <- bind_rows(lapply(split(pi, pi$group), function(dat) {
  pairwise_wilcox_table(dat$avg_pi, droplevels(dat$pop))
}), .id = "group")
write.table(pi_tests, "../figures/Nucleotide_diversity_wilcoxon.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

###################
## Heterozygosity #
###################

hetstat <- merge(hetstat,bffsff_pop, by.x = "Sample", by.y = "V1", all.x = TRUE)

h1 <- ggplot(hetstat,aes(reorder(Pop,Het_frac),Het_frac,fill="grey")) +
  geom_boxplot(outliers=FALSE,fill="grey",width=0.75)+
  theme_bw()+
  xlab("")+
  ylim(0.002,0.006)+
  theme(text = element_text(size=18),axis.text.x = element_text(angle = 45, vjust = 0.5))
hetstat$V2 <- factor(hetstat$V2, levels = c("INDO_outgroup","BFF","SFF","LRFF","GHFF", "PNG","Wet_tropics","INDO","N_AUS","NQ","E_COAST"))
h2 <- ggplot(hetstat[!is.na(hetstat$V2) & hetstat$V2!="INDO_outgroup",],aes(V2,Het_frac,fill=Pop)) +
  scale_fill_manual(values = c("lightblue","firebrick"))+
  geom_boxplot(outliers=FALSE,width=0.75)+
  xlab("")+
  ylim(0.002,0.005)+
  theme_bw()+
  theme(text = element_text(size=18),axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "none")

plot_grid(h1,h2)
ggsave("../figures/Heterozygosity_species_populations.pdf",height = 4,width = 8)

het_summary <- bind_rows(
  hetstat %>% group_by(group = Pop) %>% summarize(count = n(), mean = mean(Het_frac, na.rm = TRUE), sd = sd(Het_frac, na.rm = TRUE)),
  hetstat %>% filter(!is.na(V2) & V2 != "INDO_outgroup") %>% group_by(group = as.character(V2)) %>%
    summarize(count = n(), mean = mean(Het_frac, na.rm = TRUE), sd = sd(Het_frac, na.rm = TRUE))
)
write.table(het_summary, "../figures/Heterozygosity_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

het_reg <- hetstat[!is.na(hetstat$V2) & hetstat$V2!="INDO_outgroup",]
het_tests <- bind_rows(
  species = pairwise_wilcox_table(hetstat$Het_frac, hetstat$Pop),
  population = pairwise_wilcox_table(het_reg$Het_frac, droplevels(het_reg$V2)),
  .id = "group")
write.table(het_tests, "../figures/Heterozygosity_wilcoxon.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

##############
## Fst Plot  #
##############

# Plot FST for species
rownames(fst_species) <- fst_species$X
fst_species <- fst_species %>% select(-X)
fst_species <- as.matrix(fst_species)
fst_species<- melt(fst_species,na.rm = TRUE)

ggplot(fst_species, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                       midpoint = 0, limit = c(0,0.5), space = "Lab",
                       name="Fst") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()+
  xlab("")+
  ylab("")
ggsave("../figures/Fst_species.pdf")

# Plot FST for P. alecto and P. conspicillatus regional populations (without the
# P. alecto alecto outgroup INDO_outgroup)
rownames(fst_bffsff) <- fst_bffsff$X
fst_bffsff <- fst_bffsff %>% select(-X)
fst_bffsff <- fst_bffsff[rownames(fst_bffsff) != "INDO_outgroup", colnames(fst_bffsff) != "INDO_outgroup"]
fst_bffsff <- as.matrix(fst_bffsff)
fst_bffsff<- melt(fst_bffsff,na.rm = TRUE)

ggplot(fst_bffsff, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                       midpoint = 0, limit = c(0,0.15), space = "Lab",
                       name="Fst") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()+
  xlab("")+
  ylab("")
ggsave("../figures/Fst_bff_sff_regional.pdf")
