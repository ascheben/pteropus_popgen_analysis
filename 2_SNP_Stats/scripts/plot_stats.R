library(ggplot2)
library(dplyr)
library(cowplot)


get_box_stats <- function(y, upper_limit = 0.02) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 4), "\n",
      "Median =", round(median(y), 4), "\n"
    )
  ))
}

hetstat <- read.delim("../data/pteropus_heterozygosity.txt",sep = '\t',header=TRUE)
bffsff_pop <- read.delim("../../0_Metadata/data/palecto_pconspicillatus.popmap",sep = '\t',header=FALSE)
# Read Fst tables
fst_species <-  read.delim("../data/pteropus.fst_summary.tsv")
fst_bffsff <-  read.delim("../data/palecto_pconspicillatus.fst_summary.tsv")

# Read tab-separated popmap ("sample\tgroup") matching names in VCF

pi <- read.delim("../../2_SNP_Stats/data/pteropus_nucleotide_diversity.txt.gz",header = TRUE, sep = "\t")
pi <- pi %>% mutate(group=ifelse(pop == "BFF" | pop == "GHFF" | pop == "LRFF" | pop == "INDO_outgroup" | pop == "SFF","species","population")) %>%
  mutate(species=ifelse(pop=="PNG" | pop=="Wet_tropics","SFF","BFF"))
pi <- pi %>% filter(no_sites>=50)

pi$pop <- as.factor(pi$pop)
ggplot(pi[pi$group=="species",],aes(pop,avg_pi,fill="darkred")) + geom_boxplot(outliers = FALSE) + theme_bw()+
  theme(text = element_text(size=18))+
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9)+
  stat_summary(fun = "mean", geom = "point", shape = 2, size = 6, color = "black")

pi$pop <- factor(pi$pop, levels = c("INDO_outgroup","BFF","SFF","LRFF","GHFF", "PNG","Wet_tropics","INDO","N_AUS","NQ","E_COAST"))
p1 <- ggplot(pi[pi$group=="population",],aes(pop,avg_pi,fill=species)) + geom_boxplot(outliers = FALSE,width=0.75) + theme_bw()+
  scale_fill_manual(values = c("lightblue","firebrick"))+
  xlab("")+
  ylab("\u03c0")+
  theme(text = element_text(size=18,angle = 45),axis.text.x = element_text(angle = 45, vjust = 0.5))+
  #theme(text = element_text(size=18),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 2)+
  stat_summary(fun = "mean", geom = "point", shape = 2, size = 2, color = "black")
#ggsave("pi_100kbp_populations.pdf",height = 6,width = 6)


p2 <- ggplot(pi[pi$group=="species",],aes(reorder(pop,avg_pi),avg_pi,fill="grey")) +
  geom_boxplot(outliers = FALSE,width=0.61) + theme_bw()+
  scale_fill_manual(values = c("grey"))+
  xlab("")+
  ylab("\u03c0")+
  #theme(text = element_text(size=18),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(text = element_text(size=18),axis.text.x = element_text(angle = 45, vjust = 0.5))+
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 2)+
  stat_summary(fun = "mean", geom = "point", shape = 2, size = 2, color = "black")
#ggsave("pi_100kbp_species.pdf",height = 6,width = 5)

plot_grid(p1,p2)
ggsave("../figures/Nucleotide_diversity_all_species.pdf",height = 6,width = 11)


## Heterozygosity statistics

hetstat <- merge(hetstat,bffsff_pop, by.x = "sample", by.y = "V1", all.x = TRUE)

h1 <- ggplot(hetstat,aes(reorder(pop,het_frac),het_frac,fill="grey")) + 
  geom_boxplot(outliers=FALSE,fill="grey",width=0.75)+
  #stat_summary(fun = "mean", geom = "point", shape = 2, size =2, color = "black")+
  theme_bw()+
  xlab("")+
  ylim(0.003,0.007)+
  theme(text = element_text(size=18),axis.text.x = element_text(angle = 45, vjust = 0.5))
hetstat$V2 <- factor(hetstat$V2, levels = c("INDO_outgroup","BFF","SFF","LRFF","GHFF", "PNG","Wet_tropics","INDO","N_AUS","NQ","E_COAST"))
h2 <- ggplot(hetstat[!is.na(hetstat$V2) & hetstat$V2!="INDO_outgroup",],aes(V2,het_frac,fill=pop)) +
  scale_fill_manual(values = c("lightblue","firebrick"))+
  geom_boxplot(outliers=FALSE,width=0.75)+
  xlab("")+
  ylim(0.003,0.007)+
  #stat_summary(fun = "mean", geom = "point", shape = 2, size = 2, color = "black")+
  theme_bw()+
  theme(text = element_text(size=18),axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "none")

plot_grid(h1,h2)
ggsave("../figures/Heterozygosity_species_populations.pdf",height = 4,width = 8)

hetstat %>% group_by(pop) %>% summarize(
  count = n(),
  mean = mean(het_frac, na.rm = TRUE), 
  sd = sd(het_frac, na.rm = TRUE)) 
  
  
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
                       midpoint = 0, limit = c(0,0.4), space = "Lab",
                       name="Fst") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()+
  xlab("")+
  ylab("")
ggsave("../figures/Fst_species.pdf")

# Plot FST for BFFSFF with INDO outgroup

rownames(fst_bffsff) <- fst_bffsff$X
rownames(fst_bffsff)

fst_bffsff <- fst_bffsff %>% select(-X)
fst_bffsff <- as.matrix(fst_bffsff)
fst_bffsff<- melt(fst_bffsff,na.rm = TRUE)

ggplot(fst_bffsff, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                       midpoint = 0, limit = c(0,0.3), space = "Lab",
                       name="Fst") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()+
  xlab("")+
  ylab("")
  # Plot FST for BFFSFF without INDO outgroup
ggsave("../figures/Fst_bff_sff_regional.pdf")
