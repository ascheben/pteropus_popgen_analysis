# PCA, sampling maps, phylogenies (RAxML and CASTER-site) and isolation-by-distance
# analyses for the four flying-fox species and for P. alecto and P. conspicillatus.
#
# Run from this directory:  cd 3_SNP_Structure/scripts && Rscript pca_and_phylo_plot.R
#
# Inputs:
#   ../../1_SNP_Calling/data/pteropus.vcf.gz                        "all" SNP set
#   ../../1_SNP_Calling/data/palecto_pconspicillatus_narrow.vcf.gz  "alecto ex alecto alecto + conspicillatus" SNP set
#   ../../2_SNP_Stats/data/palecto_pconspicillatus.kinship.tsv      PLINK2 KING kinship ("alecto + conspicillatus" set)
#   ../../0_Metadata/data/{pteropus.popmap,pteropus_metadata.txt,palecto_pconspicillatus.popmap}
#   ../data/raxml/RAxML_bipartitions.{pteropus,palecto_pconspicillatus}   RAxML trees
#   ../data/raxml/CASTER.{pteropus,palecto_pconspicillatus}               CASTER-site trees (rooted with
#                                                                         P. scapulatus and P. alecto alecto)
# The VCFs are converted once to SNPRelate GDS files (../data/*.gds, not tracked by git).
#
# Outputs (../figures/):
#   Geography_lrff_ghff.pdf                          Figure 3A
#   Geography_bff_sff.pdf                            Figure 3B
#   PCA_all_species.pdf                              Figure 3C
#   PCA_PC1_PC2_bff_sff_species_regions.pdf          Figure 3D
#   Distance_from_Sawu_by_PC1.pdf                    Figure 4C
#   GeneticDistance_by_GeoDistance.pdf               Figure 4D
#   Geography_bff_sff_regions.pdf                    Figure S1
#   PCA_PC1_PC2_bff_sff_species_regions_nomiss.pdf   Figure S9 (SNPs without missing genotypes)
#   Phylogram_bff_sff.pdf                            Figure S12 (RAxML)
#   Phylogram_lrff_ghff.pdf                          Figure S13 (RAxML)
#   Phylogram_caster_site_bff_sff.pdf                Figure S14 (CASTER-site)
#   Phylogram_caster_site_lrff_ghff.pdf              Figure S15 (CASTER-site)
#   Species_PCA.tsv, PCA_bff_sff.tsv                 PC coordinates (Table S11)
#   mantel_test_genetic_vs_geographic_distance.txt   Mantel test, P. alecto (Table S17)
#   BFF_genetic_distance.tsv                         pairwise genetic and geographic distances, P. alecto
#   PCA_PC1_PC3_*, Kinship_by_distance.pdf, Geography_all_species.pdf, Phylogram_all_species.pdf
#                                                    supporting plots (not in manuscript)
# Labels, colours of tip labels, legends and panel layout were finalised in a graphics editor.

library(SNPRelate)
library(gdsfmt)
library(ggrepel)
library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(colorspace)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(geosphere)
library(adegenet)
library(vegan)
library(vcfR)
library(poppr)
library(tidyr)
library(tibble)
library(ggpointdensity)
library(viridis)

ff_vcf <- "../../1_SNP_Calling/data/pteropus.vcf.gz"
bffsff_vcf <- "../../1_SNP_Calling/data/palecto_pconspicillatus_narrow.vcf.gz"
bffsff_kinship <- read.delim("../../2_SNP_Stats/data/palecto_pconspicillatus.kinship.tsv",sep='\t')

# Read tab-separated popmap ("sample\tgroup") matching names in VCF
popfile <- read.csv("../../0_Metadata/data/pteropus.popmap", header = FALSE, sep = "\t")
geography <- read.csv("../../0_Metadata/data/pteropus_metadata.txt", header = TRUE, sep = "\t")
bffsff_popfile <- read.csv("../../0_Metadata/data/palecto_pconspicillatus.popmap", header = FALSE, sep = "\t")

# Store tree locations
species_raxml_tree <- "../data/raxml/RAxML_bipartitions.pteropus"
bffsff_raxml_tree <- "../data/raxml/RAxML_bipartitions.palecto_pconspicillatus"
species_caster_tree <- "../data/raxml/CASTER.pteropus"
bffsff_caster_tree <- "../data/raxml/CASTER.palecto_pconspicillatus"

dir.create("../figures", showWarnings = FALSE)

# Open a SNPRelate GDS file, converting the VCF on first use
open_gds <- function(vcf, gdsfile) {
  if (!file.exists(gdsfile)) snpgdsVCF2GDS(vcf, gdsfile, method="biallelic.only")
  snpgdsOpen(gdsfile, allow.duplicate=TRUE)
}

# Data frame of the first four PCs and axis labels with % variance explained
pca_table <- function(pca) {
  round.pc.percent <- round(pca$varprop*100, 2)
  tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2], # the second eigenvector
                    EV3 = pca$eigenvect[,3], # the third eigenvector
                    EV4 = pca$eigenvect[,4], # the fourth eigenvector
                    stringsAsFactors = FALSE)
  tab$sample.id <- sub(".*__", "", tab$sample.id)
  list(tab = tab,
       lab = paste0("PC", 1:4, " (", round.pc.percent[1:4], "%)"))
}

#####################################
## PCA of all species (Figure 3C)  ##
#####################################

ff_cols <- rainbow_hcl(5)
# Add popmap header
colnames(popfile) <- c("sample", "species","region","location")
genofile <- open_gds(ff_vcf, "../data/ff_out.gds")
# missing.rate = NaN: no missingness filter (the default in SNPRelate 1.40 used for
# the manuscript; newer versions exclude SNPs with > 1% missing genotypes by default)
pca <- pca_table(snpgdsPCA(genofile,autosome.only=FALSE,num.thread=4,missing.rate=NaN))
snpgdsClose(genofile)
tabpops <- merge(pca$tab,popfile, by.x = "sample.id" , by.y = "sample" )

genogeo <- merge(tabpops,geography,by.x = c("sample.id"),by.y = c("Sample.Identifier"))
# Indonesian P. alecto separated on PC1 are P. alecto alecto (IBFF)
genogeo <- genogeo %>% mutate(Species2 = ifelse(species == "BFF" & Country == "Indonesia" & EV1>0,"IBFF",species))

shapes <- c("BFF" = 21, "GHFF" = 24, "IBFF" = 25, "LRFF" = 23,"SFF" = 22)
ggplot(genogeo,aes(x=EV1,y=EV2,fill=Species2,shape=Species2,label=sample.id))+
  geom_point(size=3,alpha=0.5) +
  scale_shape_manual(values=shapes)+
  ylim(-0.1,0.1)+
  xlim(-0.1,0.1)+
  theme_bw()+
  theme(text = element_text(size=18)) +
  scale_color_manual(values=c(ff_cols)) +
  xlab(pca$lab[1]) +
  ylab(pca$lab[2])
ggsave("../figures/PCA_all_species.pdf",height=7,width=8)

genogeo_tab <- genogeo %>% select(sample.id,EV1,EV2,EV3,EV4)
write.table(genogeo_tab,"../figures/Species_PCA.tsv",quote = F,row.names = F,sep = '\t')

###############################################
## Sampling maps (Figures 3A, 3B and S1)     ##
###############################################

gg_uniq <- genogeo %>% group_by(Latitude,Longitude,Species2)  %>% distinct(pick(Latitude,Longitude,Species2))
gg_uniq_regions <- merge(genogeo,bffsff_popfile, by.x = c("sample.id"),by.y = c("V1")) %>% group_by(Latitude,Longitude,V2)  %>% distinct(pick(Latitude,Longitude,V2,Species2))

ff_map <- ne_states(country = c("indonesia","australia","papua new guinea"), returnclass = "sf")
ggplot() +
  geom_sf(data = ff_map,color = "grey80")+
  geom_point(data = gg_uniq[gg_uniq$Species2=="BFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 21, fill = "#E495A5",alpha=1)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="IBFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 25, fill = "#65BC8C",alpha=1)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="SFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 22, fill = "#C29DDE",alpha=1)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="LRFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 23, fill = "#55B8D0",alpha=1)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="GHFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 24, fill = "#BDAB66",alpha=1)+
  coord_sf(xlim = c(95, 160), ylim = c(5, -45), expand = FALSE)+
  theme_bw()
ggsave("../figures/Geography_all_species.pdf")

ggplot() +
  geom_sf(data = ff_map,color = "grey80")+
  geom_point(data = gg_uniq[gg_uniq$Species2=="BFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 21, fill = "#E495A5",alpha=0.5)+
  geom_jitter(data = gg_uniq[gg_uniq$Species2=="IBFF",], aes(x = Longitude, y = Latitude), size = 2,
              shape = 25, fill = "#65BC8C",alpha=0.5)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="SFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 22, fill = "#C29DDE",alpha=0.5)+
  coord_sf(xlim = c(95, 160), ylim = c(5, -45), expand = FALSE)+
  theme_bw()
ggsave("../figures/Geography_bff_sff.pdf")

ggplot() +
  geom_sf(data = ff_map,color = "grey80")+
  geom_point(data = gg_uniq[gg_uniq$Species2=="LRFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 23, fill = "#55B8D0",alpha=0.5)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="GHFF",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 24, fill = "#BDAB66",alpha=0.5)+
  coord_sf(xlim = c(95, 160), ylim = c(5, -45), expand = FALSE)+
  theme_bw()
ggsave("../figures/Geography_lrff_ghff.pdf")

ggplot() +
  geom_sf(data = ff_map,color = "grey80")+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="E_COAST",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 21, fill = "#65BC8C",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="INDO",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 21, fill = "#beaed4",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="N_AUS",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 21, fill = "#ffff99",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="NQ",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 21, fill = "#386cb0",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="PNG",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 22, fill = "#f0027f",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="Wet_tropics",], aes(x = Longitude, y = Latitude), size = 2,
             shape = 22, fill = "#bf5b17",alpha=0.5)+
  coord_sf(xlim = c(95, 160), ylim = c(5, -45), expand = FALSE)+
  theme_bw()
ggsave("../figures/Geography_bff_sff_regions.pdf")

##################################################################
## PCA of P. alecto and P. conspicillatus (Figures 3D and S9)   ##
##################################################################

# Use only BFF and SFF
colnames(bffsff_popfile) <- c("sample", "group")
genofile <- open_gds(bffsff_vcf, "../data/out_bff_sff_narrow.gds")

# Figure S9: only SNPs without missing genotypes (216 SNPs)
snp_missing <- snpgdsSNPRateFreq(genofile, with.snp.id = TRUE)
nomiss_snps <- snp_missing$snp.id[snp_missing$MissingRate == 0]

bffsff_pca <- function(pca) {
  genogeo_bffsff <- merge(pca$tab,bffsff_popfile, by.x = "sample.id" , by.y = "sample" )
  merge(genogeo_bffsff,geography,by.x = c("sample.id"), by.y = c("Sample.Identifier") )
}

plot_bffsff_pca <- function(genogeo_bffsff, lab, outfile) {
  ggplot(genogeo_bffsff,aes(x=EV1,y=EV2,fill=BFF_SFF_Population_ID,shape=Species, label=sample.id))+
    scale_shape_manual(values = c(21,22))+
    scale_fill_manual(values=c(c("#65BC8C","#beaed4","#ffff99","#386cb0","#f0027f","#bf5b17","#fdc086"))) +
    geom_point(size=4,alpha=0.5) +
    guides(fill = guide_legend(override.aes = list(shape=21)))+
    theme_bw()+
    theme(text = element_text(size=18)) +
    xlab(lab[1]) +
    ylab(lab[2])
  ggsave(outfile,height=3.5,width=7.5)
}

pca_nomiss <- pca_table(snpgdsPCA(genofile,autosome.only=FALSE,snp.id=nomiss_snps,num.thread=2,missing.rate=NaN))
plot_bffsff_pca(bffsff_pca(pca_nomiss), pca_nomiss$lab, "../figures/PCA_PC1_PC2_bff_sff_species_regions_nomiss.pdf")

# Figure 3D: all SNPs
pca <- pca_table(snpgdsPCA(genofile,autosome.only=FALSE,num.thread=2,missing.rate=NaN))
snpgdsClose(genofile)
genogeo_bffsff <- bffsff_pca(pca)
plot_bffsff_pca(genogeo_bffsff, pca$lab, "../figures/PCA_PC1_PC2_bff_sff_species_regions.pdf")

ggplot(genogeo_bffsff,aes(x=EV1,y=EV3,fill=BFF_SFF_Population_ID,shape=Species, label=sample.id))+
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values=c(c("#65BC8C","#beaed4","#ffff99","#386cb0","#f0027f","#bf5b17","#fdc086"))) +
  geom_point(size=4,alpha=0.5) +
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  theme_bw()+
  theme(text = element_text(size=18)) +
  xlab(pca$lab[1]) +
  ylab(pca$lab[3])
ggsave("../figures/PCA_PC1_PC3_bff_sff_species_regions.pdf",height=3.5,width=7.5)

ggplot(genogeo_bffsff,aes(x=EV1,y=EV3,fill=BFF_SFF_Population_ID,shape=Species, label=sample.id))+
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values=c(c("#65BC8C","#beaed4","#ffff99","#386cb0","#f0027f","#bf5b17","#fdc086"))) +
  geom_point(size=4,alpha=0.5) +
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  geom_text_repel(size=2) +
  theme_bw()+
  theme(text = element_text(size=18)) +
  xlab(pca$lab[1]) +
  ylab(pca$lab[3])
ggsave("../figures/PCA_PC1_PC3_bff_sff_species_regions_label.pdf",height=7,width=12)

genogeo_bffsff_tab <- genogeo_bffsff %>% select(sample.id,EV1,EV2,EV3,EV4)
write.table(genogeo_bffsff_tab,"../figures/PCA_bff_sff.tsv",quote = F,row.names = F,sep = '\t')

##############################################################
## Phylogenies (RAxML: Figures S12, S13; CASTER: S14, S15)  ##
##############################################################

specieslist <- split(genogeo$sample.id, genogeo$Species2)

# Region of the roost for colouring the P. poliocephalus and P. scapulatus trees
loc2regiondf <- data.frame(Roost = c("Canungra","Adelaide","Katherine","Tolga","Isa","Wingham","Pt Macquarie","Casino","MacLean","Townsville","Regents Park","Charters","Broome","Lismore","Kyogle","Bellingen","Bundgeam","Ingham","Pt Douglas"),
                           custom_region = c("ECOAST","SAUS","N_AUS","Wet_tropics","NQ","ECOAST","ECOAST","ECOAST","ECOAST","NQ","ECOAST","NQ","N_AUS","ECOAST","ECOAST","ECOAST","ECOAST","NQ","Wet_tropics"))
genogeo$Roost <- trimws(genogeo$Roost)
genogeo <- merge(genogeo,loc2regiondf,by = "Roost", all.x =TRUE)
regionlist <- split(genogeo$sample.id, genogeo$custom_region)

bffsff_list <- split(bffsff_popfile$sample, bffsff_popfile$group)

# Species tree coloured by region with the P. alecto / P. conspicillatus clades collapsed
plot_lrff_ghff_tree <- function(tree, collapse_nodes, outfile) {
  p <- ggtree(groupOTU(tree, regionlist), aes(color=group), layout="rectangular",branch.length="none") + geom_tiplab(size=1) +
    scale_color_manual(values=c(rainbow_hcl(7))) +
    theme(legend.position="right")+
    geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label))& as.numeric(label) > 70))+
    scale_size_manual(values=c(1, .1))
  for (n in collapse_nodes) p <- ggtree::collapse(p, node=n)
  print(p)
  ggsave(outfile,height=10,width=8)
}

# RAxML species tree, rooted with P. scapulatus (species assignment from the popmap,
# which accounts for re-determined samples such as 1510019_LRFF)
sp_raxml <- read.tree(species_raxml_tree)
lrff_samplesids <- intersect(popfile$sample[popfile$species == "LRFF"], sp_raxml$tip.label)
sp_raxml <- root(unroot(sp_raxml), outgroup = lrff_samplesids, resolve.root = TRUE)

ggtree(groupOTU(sp_raxml, specieslist), aes(color=group),branch.length="none", layout="rectangular") +
  scale_color_manual(values=c(rainbow_hcl(6))) +
  theme(legend.position="right")+
  geom_tiplab(size =0.5)+
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1))  +
  ggplot2::xlim(0, 60)
ggsave("../figures/Phylogram_all_species.pdf",height=10,width=8)

plot_lrff_ghff_tree(sp_raxml,
                    getMRCA(sp_raxml, popfile$sample[popfile$species %in% c("BFF","SFF")]),
                    "../figures/Phylogram_lrff_ghff.pdf")

# RAxML P. alecto and P. conspicillatus tree coloured by region
ggtree(groupOTU(read.tree(bffsff_raxml_tree), bffsff_list), aes(color=group),branch.length="none", layout="rectangular") +
  scale_color_manual(values=c("grey","#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17")) +
  theme(legend.position="right")+
  geom_tiplab(size =0.75)+
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1))  +
  ggplot2::xlim(0, 60)
ggsave("../figures/Phylogram_bff_sff.pdf",height=10,width=8)

# CASTER-site species tree (rooted with P. scapulatus): collapse P. alecto alecto and
# the P. alecto + P. conspicillatus clade
sp_caster <- read.tree(species_caster_tree)
caster_species <- geography$Species[match(sp_caster$tip.label, geography$Sample.Identifier)]
plot_lrff_ghff_tree(sp_caster,
                    c(getMRCA(sp_caster, sp_caster$tip.label[caster_species == "IBFF"]),
                      getMRCA(sp_caster, sp_caster$tip.label[caster_species %in% c("BFF","SFF")])),
                    "../figures/Phylogram_caster_site_lrff_ghff.pdf")

# CASTER-site P. alecto and P. conspicillatus tree (rooted with P. alecto alecto), with the
# two daughter clades of the P. alecto + P. conspicillatus ancestor flipped (P. conspicillatus on top)
bffsff_caster <- read.tree(bffsff_caster_tree)
caster_species <- geography$Species[match(bffsff_caster$tip.label, geography$Sample.Identifier)]
bffsff_node <- getMRCA(bffsff_caster, bffsff_caster$tip.label[caster_species %in% c("BFF","SFF")])
flip_nodes <- bffsff_caster$edge[bffsff_caster$edge[,1] == bffsff_node, 2]
flip(ggtree(groupOTU(bffsff_caster, bffsff_list), aes(color=group),branch.length="none", layout="rectangular") +
  scale_color_manual(values=c("grey","#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17")) +
  theme(legend.position="right")+
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1))  +
  ggplot2::xlim(0, 65)+
  geom_tiplab(size =0.75), flip_nodes[2], flip_nodes[1])
ggsave("../figures/Phylogram_caster_site_bff_sff.pdf",height=10,width=8)

####################################################################
## Isolation by distance in P. alecto (Figures 4C, 4D; Table S17) ##
####################################################################

genogeo_bffsff <- genogeo_bffsff %>% rowwise() %>% mutate(km_distance_from_sawu=distm(c(121.9167,-10.4833), c(Longitude,Latitude), fun = distHaversine)/1000)

ggplot(genogeo_bffsff,aes(km_distance_from_sawu,EV1,fill=Species,shape=Species)) +
  geom_point(size=4,alpha=0.5)+
  scale_fill_manual(values=c("#E495A5","#C29DDE"))+
  scale_shape_manual(values=c(21,22))+
  ylab("PC1")+
  xlab("Geographic distance from Sawu (km)")+
  theme_bw()+
  theme(text = element_text(size=18),legend.title = element_blank())
ggsave("../figures/Distance_from_Sawu_by_PC1.pdf",width =7,height = 4)

# Pairs of P. alecto samples with kinship and geographic distance
distance_geo_gen <- merge(bffsff_kinship,geography, by.x='IID1', by.y='Sample.Identifier')
distance_geo_gen <- merge(distance_geo_gen,geography, by.x='IID2', by.y='Sample.Identifier')
distance_geo_gen_bff <- distance_geo_gen %>%
  select(IID1,IID2,KINSHIP,Species.x,Species.y,Longitude.x,Latitude.x,Longitude.y,Latitude.y,
         BFF_SFF_Population_ID.x,BFF_SFF_Population_ID.y) %>%
  filter(Species.y == 'BFF',Species.x == 'BFF') %>%
  rowwise() %>%
  mutate(km_distance=distm(c(Longitude.x,Latitude.x), c(Longitude.y,Latitude.y), fun = distHaversine)/1000) %>%
  ungroup() %>%
  mutate(indo=ifelse((BFF_SFF_Population_ID.x == "INDO" & BFF_SFF_Population_ID.y !="INDO")| (BFF_SFF_Population_ID.y == "INDO" & BFF_SFF_Population_ID.x !="INDO"),"Indonesia-Australia","Australia-Australia"))

ggplot(distance_geo_gen_bff,aes(km_distance,KINSHIP,color=indo,fill=indo)) + geom_point(shape=21,size=2,alpha=0.5)+
  geom_smooth(method='lm', formula= y~x,color="black")+
  scale_fill_manual(values=c("darkred","darkblue"))+
  scale_color_manual(values=c("darkred","darkblue"))+
  ylab("Kinship coefficient")+
  xlab("Geographic distance (km)")+
  theme_bw()+
  theme(text = element_text(size=18),legend.title = element_blank())
ggsave("../figures/Kinship_by_distance.pdf",width =7,height = 4)

## Mantel test of genetic vs geographic distance (86 P. alecto samples)

mantel_samples <- unique(c(distance_geo_gen_bff$IID1,distance_geo_gen_bff$IID2))

# Genetic distance (proportion of allelic differences, poppr::diss.dist) from the
# SNPs without missing genotypes (the same 216 SNPs as Figure S9)
dist_vcf <- read.vcfR(bffsff_vcf, verbose = FALSE)
dist_vcf <- dist_vcf[rowSums(is.na(extract.gt(dist_vcf))) == 0, ]
gl_data <- vcfR2genind(dist_vcf[, c("FORMAT", mantel_samples)])
genetic_dist_matrix <- poppr::diss.dist(x = gl_data, percent = TRUE, mat = TRUE)

# Geographic distance matrix in the same sample order
sample_order <- rownames(genetic_dist_matrix)
geo_dist_matrix <- matrix(0, nrow = length(sample_order), ncol = length(sample_order),
                          dimnames = list(sample_order, sample_order))
for (i in 1:nrow(distance_geo_gen_bff)) {
  id1 <- distance_geo_gen_bff$IID1[i]
  id2 <- distance_geo_gen_bff$IID2[i]
  geo_dist_matrix[id1, id2] <- distance_geo_gen_bff$km_distance[i]
  geo_dist_matrix[id2, id1] <- distance_geo_gen_bff$km_distance[i]
}

set.seed(1)
IBD <- vegan::mantel(as.dist(genetic_dist_matrix), as.dist(geo_dist_matrix), method="pearson")
IBD

writeLines(c(
  "--- Mantel Test Results (Isolation by Distance) ---",
  paste("Method:", IBD$method),
  paste("Correlation (r) statistic:", round(IBD$statistic, 4)),
  paste("Significance (p-value):", format.pval(IBD$signif, digits = 4)),
  paste("Permutations:", IBD$permutations),
  "---------------------------------------------------"
), "../figures/mantel_test_genetic_vs_geographic_distance.txt")

# Pairwise genetic vs geographic distance (Figure 4D)
genetic_dist_long <- as.data.frame(genetic_dist_matrix) %>%
  rownames_to_column("IID1") %>%
  pivot_longer(!IID1 ,names_to = "IID2", values_to = "GenDist")
distance_geo_gen_bff_pop <- merge(distance_geo_gen_bff, genetic_dist_long, by = c('IID1','IID2'))

ggplot(distance_geo_gen_bff_pop,aes(km_distance,GenDist)) +
  geom_pointdensity(size = 1.2)+
  scale_color_viridis()+
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    color = "black",
    fill = "grey70", # Lighter fill for SE band
    alpha = 0.3,
    linewidth = 1 # Make the line slightly thicker
  ) +
  ylab("Genetic distance")+
  xlab("Geographic distance (km)")+
  theme_bw()+
  theme(text = element_text(size=18),legend.title = element_blank())
ggsave("../figures/GeneticDistance_by_GeoDistance.pdf")

write.table(distance_geo_gen_bff_pop,"../figures/BFF_genetic_distance.tsv",sep='\t',quote = F,row.names = F)
