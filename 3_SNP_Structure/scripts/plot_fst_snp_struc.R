library(SNPRelate)
library(gdsfmt)
library(gridExtra)
library(ggrepel)
library(ggtree)
library(ape)
library(reshape2)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(colorspace)
library(cowplot)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)


ff_vcf <- "1_SNP_Calling\\data\\bcftools.no_5bp_flanks.biallelic.whitelist.snps_mim09_biallelic_minDP5_mm05_maf005.indepr03.vcf.gz"
bffsff_vcf <- "1_SNP_Calling\\data\\bcftools.no_5bp_flanks.biallelic.whitelist.snps_mim09_biallelic_minDP5.bffsff_mm05_maf005.indepr03.vcf.gz"

# Read tab-separated popmap ("sample\tgroup") matching names in VCF
popfile <- read.csv("0_Metadata\\popmap_all_fix.txt", header = FALSE, sep = "\t")
geography <- read.csv("0_Metadata\\bat_metadata_sheet_all.txt", header = TRUE, sep = "\t") 
bffsff_popfile <- read.csv("0_Metadata\\popmap_bffsff_big_groups.txt", header = FALSE, sep = "\t")

# Read Fst tables
fst_species <-  read.delim("2_SNP_Stats\\species\\bcftools.no_5bp_flanks.biallelic.whitelist.snps_mim09_biallelic_minDP5_mm05_maf005.indepr03.p.fst_summary.tsv")
fst_bffsff <-  read.delim("2_SNP_Stats\\BFFSFF_custom\\bcftools.no_5bp_flanks.biallelic.whitelist.snps_mim09_biallelic_minDP5.bffsff_mm05_maf005.indepr03.p.fst_summary.tsv")

# Store tree locations
#species_tree <- "3_SNP_Structure/data/raxml/RAxML_bipartitions.bcftools.no_5bp_flanks.biallelic.whitelist.snps_mim09_biallelic_minDP5_mm05_maf005.indepr03"
species_tree <- "3_SNP_Structure/data/raxml/raxml_species_rerooted"
#bffsff_tree <- "3_SNP_Structure/data/raxml/RAxML_bipartitions.bcftools.no_5bp_flanks.biallelic.whitelist.snps_mim09_biallelic_minDP5.bffsff_mm05_maf005.indepr03"
#bffsff_tree <- "3_SNP_Structure/data/raxml/raxml_bffsff_rerooted"
bffsff_tree <- "3_SNP_Structure/data/raxml/RAxML_bipartitions.bcftools.no_5bp_flanks.biallelic.whitelist.snps_mim09_biallelic_minDP5.bffsff_mm05_maf005.indepr03_outgroup"



##############
## Figure 1  #
##############
ff_cols <- rainbow_hcl(5)


# Add popmap header
colnames(popfile) <- c("sample", "species","region","location")
# Set output file names
outname <- "ff_out"
gdsfile <- paste(outname, "gds", sep=".")
pdffile1 <- paste(outname, "_label_pca.pdf", sep="")
pdffile2 <- paste(outname, "_pca.pdf", sep="")
# Convert VCF to GDS
#snpgdsVCF2GDS(ff_vcf,gdsfile, method="biallelic.only")
# Read in genotypes
genofile <- snpgdsOpen(gdsfile, allow.duplicate=TRUE)
# LD pruning option (turned off by default)
#set.seed(1000)
#snpset <- snpgdsLDpruning(genofile,autosome.only=FALSE, slide.max.bp = 1)
# PCA 
#snpset.id <- unlist(snpset)
#pca <- snpgdsPCA(genofile,autosome.only=FALSE, snp.id=snpset.id,num.thread=2)
pca <- snpgdsPCA(genofile,autosome.only=FALSE,num.thread=4)

pc.percent <- pca$varprop*100
round.pc.percent <- round(pc.percent, 2)
# Get % explained and store in vectors below
ev1pc <- round.pc.percent[1]
ev2pc <- round.pc.percent[2]
ev1pc <- paste("PC1 (", ev1pc, "%)", sep ="")
ev2pc <- paste("PC2 (", ev2pc, "%)", sep ="")
ev3pc <- round.pc.percent[3]
ev4pc <- round.pc.percent[4]
ev3pc <- paste("PC3 (", ev3pc, "%)", sep ="")
ev4pc <- paste("PC4 (", ev4pc, "%)", sep ="")
# Prepare PCA results
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1], # the first eigenvector
                  EV2 = pca$eigenvect[,2], # the second eigenvector
                  EV3 = pca$eigenvect[,3], # the third eigenvector
                  EV4 = pca$eigenvect[,4], # the fourth eigenvector
                  stringsAsFactors = FALSE)
tab$sample.id <- sub(".*__", "", tab$sample.id)
tabpops <- merge(tab,popfile, by.x = "sample.id" , by.y = "sample" )
# Set plot theme
theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))


## Analysis 16.05.2024
#geography <- read.delim("bat_metadata_geography.txt",header = T,sep='\t')
#longlat <- read.delim("bat_long_lat.txt",header=T,sep='\t')
#genogeo <- merge(tab,geography, by.x = "sample.id" , by.y = "name" )
#genogeo <- merge(genogeo,longlat, by.x = "sample.id" , by.y = "sample.id", all.x = T)


genogeo <- merge(tabpops,geography,by.x = c("sample.id"),by.y = c("Sample_name"))
genogeo <- genogeo %>% mutate(Species2 = ifelse(species == "BFF" & Country == "Indonesia" & EV1>0,"IBFF",species))

options(ggrepel.max.overlaps = 100)
shapes <- c("BFF" = 21, "GHFF" = 24, "IBFF" = 25, "LRFF" = 23,"SFF" = 22)
ggplot(genogeo,aes(x=EV1,y=EV2,fill=Species2,shape=Species2,label=sample.id))+ 
  geom_point(size=3,alpha=0.5) + 
  scale_shape_manual(values=shapes)+
  ylim(-0.1,0.1)+
  xlim(-0.1,0.1)+
  #geom_text_repel() + 
  theme_bw()+
  theme(text = element_text(size=18)) + 
  scale_color_manual(values=c(ff_cols)) + 
  xlab(ev1pc) + 
  ylab(ev2pc)
ggsave("PCA_by_species_7by8.pdf",height=7,width=8)
ggsave("PCA_by_species.pdf")

genogeo_tab <- genogeo %>% select(sample.id,EV1,EV2,EV3,EV4)
write.table(genogeo_tab,"Species_PCA.tsv",quote = F,row.names = F,sep = '\t')

##############
## Figure 2  #
##############

# Plot geography

gg_uniq <- genogeo %>% group_by(Lat,Lon,Species2)  %>% distinct(pick(Lat,Lon,Species2))
gg_uniq_regions <- merge(genogeo,bffsff_popfile, by.x = c("sample.id"),by.y = c("V1")) %>% group_by(Lat,Lon,V2)  %>% distinct(pick(Lat,Lon,V2,Species2))

ff_map <- ne_states(country = c("indonesia","australia","papua new guinea"), returnclass = "sf")
ggplot() +
  geom_sf(data = ff_map,color = "grey80")+ 
  geom_point(data = gg_uniq[gg_uniq$Species2=="BFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 21, fill = "#E495A5",alpha=1)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="IBFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 25, fill = "#65BC8C",alpha=1)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="SFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 22, fill = "#C29DDE",alpha=1)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="LRFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 23, fill = "#55B8D0",alpha=1)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="GHFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 24, fill = "#BDAB66",alpha=1)+
  coord_sf(xlim = c(95, 160), ylim = c(5, -45), expand = FALSE)+
  theme_bw()
  #theme(panel.background = element_rect(fill = "aliceblue"))
ggsave("Geography_by_species.pdf")

ggplot() +
  geom_sf(data = ff_map,color = "grey80")+ 
  geom_point(data = gg_uniq[gg_uniq$Species2=="BFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 21, fill = "#E495A5",alpha=0.5)+
  geom_jitter(data = gg_uniq[gg_uniq$Species2=="IBFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 25, fill = "#65BC8C",alpha=0.5)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="SFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 22, fill = "#C29DDE",alpha=0.5)+
  coord_sf(xlim = c(95, 160), ylim = c(5, -45), expand = FALSE)+
  theme_bw()
ggsave("Geography_bffsff.pdf")

ggplot() +
  geom_sf(data = ff_map,color = "grey80")+ 
  geom_point(data = gg_uniq[gg_uniq$Species2=="LRFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 23, fill = "#55B8D0",alpha=0.5)+
  geom_point(data = gg_uniq[gg_uniq$Species2=="GHFF",], aes(x = Lon, y = Lat), size = 2, 
             shape = 24, fill = "#BDAB66",alpha=0.5)+
  coord_sf(xlim = c(95, 160), ylim = c(5, -45), expand = FALSE)+
  theme_bw()
ggsave("Geography_lrff_ghff.pdf")

#bffsff_cols <- c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17")

ggplot() +
  geom_sf(data = ff_map,color = "grey80")+ 
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="E_COAST",], aes(x = Lon, y = Lat), size = 2, 
             shape = 21, fill = "#65BC8C",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="INDO",], aes(x = Lon, y = Lat), size = 2, 
             shape = 21, fill = "#beaed4",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="N_AUS",], aes(x = Lon, y = Lat), size = 2, 
             shape = 21, fill = "#ffff99",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="NQ",], aes(x = Lon, y = Lat), size = 2, 
             shape = 21, fill = "#386cb0",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="PNG",], aes(x = Lon, y = Lat), size = 2, 
             shape = 22, fill = "#f0027f",alpha=0.5)+
  geom_point(data = gg_uniq_regions[gg_uniq_regions$V2=="Wet_tropics",], aes(x = Lon, y = Lat), size = 2, 
             shape = 22, fill = "#bf5b17",alpha=0.5)+
  coord_sf(xlim = c(95, 160), ylim = c(5, -45), expand = FALSE)+
  theme_bw()
ggsave("Geography_bffsff_regions.pdf")

##############
## Figure 3  #
##############

# Use only BFF and SFF

bffsff_pop <- read.csv("0_Metadata\\popmap_bffsff_big_groups.txt", header = FALSE, sep = "\t")
# Add popmap header
colnames(bffsff_pop) <- c("sample", "group")
# Set output file names
outname <- "out_bff_sff"
gdsfile <- paste(outname, "gds", sep=".")
# Convert VCF to GDS
#snpgdsVCF2GDS(bffsff_vcf,gdsfile, method="biallelic.only")
# Read in genotypes
genofile <- snpgdsOpen(gdsfile, allow.duplicate=TRUE)

# LD pruning option (turned off by default)
#set.seed(1000)
snpset <- snpgdsLDpruning(genofile,autosome.only=FALSE, slide.max.bp = 1)
#snpset <- snpgdsLDpruning(genofile,autosome.only=FALSE)

# PCA 
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile,autosome.only=FALSE, snp.id=snpset.id,num.thread=2)
pc.percent <- pca$varprop*100
round.pc.percent <- round(pc.percent, 2)
# Get % explained and store in vectors below
ev1pc <- round.pc.percent[1]
ev2pc <- round.pc.percent[2]
ev1pc <- paste("PC1 (", ev1pc, "%)", sep ="")
ev2pc <- paste("PC2 (", ev2pc, "%)", sep ="")
ev3pc <- round.pc.percent[3]
ev4pc <- round.pc.percent[4]
ev3pc <- paste("PC3 (", ev3pc, "%)", sep ="")
ev4pc <- paste("PC4 (", ev4pc, "%)", sep ="")
# Prepare PCA results
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1], # the first eigenvector
                  EV2 = pca$eigenvect[,2], # the second eigenvector
                  EV3 = pca$eigenvect[,3], # the third eigenvector
                  EV4 = pca$eigenvect[,4], # the fourth eigenvector
                  stringsAsFactors = FALSE)
tab$sample.id <- sub(".*__", "", tab$sample.id)
tabpops_bffsff <- merge(tab,bffsff_pop, by.x = "sample.id" , by.y = "sample" )

genogeo_bffsff <- merge(tabpops_bffsff,geography,by.x = c("sample.id"), by.y = c("Sample_name") )

#############################################################################################################
### WARNING : EV2 THRESHOLD MUST BE MANUALLY CHECKED USING PLOT BELOW TO ENSURE IT ONLY DETECTS TRUE IBFF ###
#############################################################################################################

genogeo_bffsff <- genogeo_bffsff %>% mutate(Species2 = ifelse(Species == "BFF" & Country == "Indonesia" & EV2< -0.1,"IBFF",Species))


genogeo_bffsff <- merge(genogeo_bffsff,bffsff_popfile, by.x = c("sample.id"),by.y = c("V1"))

  
#bffsff_cols <- c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17")
max.overlaps <- 100
ggplot(genogeo_bffsff,aes(x=EV1,y=EV3,fill=V2,shape=Species2, label=sample.id))+ 
  scale_shape_manual(values = c(21,25,22))+ 
  scale_fill_manual(values=c(c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17"))) + 
  geom_point(size=4,alpha=0.5) + 
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  #geom_text_repel() + 
  theme_bw()+
  theme(text = element_text(size=18)) + 
  xlab(ev1pc) + 
  ylab(ev3pc)
ggsave("pca_bff_sff_species_country.pdf")

ggplot(genogeo_bffsff,aes(x=EV1,y=EV2,fill=V2,shape=Species2, label=sample.id))+ 
  scale_shape_manual(values = c(21,25,22))+ 
  scale_fill_manual(values=c(c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17"))) + 
  geom_point(size=4,alpha=0.5) + 
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  #geom_text_repel() + 
  theme_bw()+
  theme(text = element_text(size=18)) + 
  xlab(ev1pc) + 
  ylab(ev2pc)

genogeo_bffsff_tab <- genogeo_bffsff %>% select(sample.id,EV1,EV2,EV3,EV4)
write.table(genogeo_bffsff_tab,"BFFSFF_PCA.tsv",quote = F,row.names = F,sep = '\t')

ggplot(genogeo_bffsff,aes(x=EV1,y=EV3,fill=V2,shape=Species2, label=sample.id))+ 
  scale_shape_manual(values = c(21,25,22))+ 
  scale_fill_manual(values=c(c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17"))) + 
  geom_point(size=4,alpha=0.5) + 
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  geom_text_repel(size=2) + 
  theme_bw()+
  theme(text = element_text(size=18)) + 
  xlab(ev1pc) + 
  ylab(ev3pc)
ggsave("pca_bff_sff_species_country_label.pdf")

##############
## Figure 4  #
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
ggsave("Fst_species.pdf")

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
ggplot(fst_bffsff[fst_bffsff$Var1!="INDO_outgroup" & fst_bffsff$Var2!="INDO_outgroup",], aes(Var2, Var1, fill = value))+
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
ggsave("Fst_bffsff_regional.pdf")

# Plot fst for finer regions

#fst <- read.delim("C:/Users/Armin/ws/projects/bat23/stats/bff_sff_region_fst.tsv")
#rownames(fst) <- fst$X
#region_names <- c(fst$X,c("Madang"))
#fst <- fst[region_names]
#fst <- as.matrix(fst)
#fst <- melt(fst,na.rm = TRUE)
#ggplot(fst, aes(x=Var1, y=Var2, fill=value)) + 
#  geom_tile()

#region_order <- c("Wet_tropics","Qld","NSW","NT","Manton_Dam","WA","Kununurra","Bigge_Island","Pt_Moresby","Madang","Sumba","Rinca","Sulawesi","Sawu")
#require(dplyr)
#fst %>%
##  mutate(Var1 = factor(Var1, levels=region_order),
#         Var2 = factor(Var2, levels=region_order)) %>%
#  ggplot(aes(Var2, Var1, fill = value))+
#  geom_tile(color = "white")+
#  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
#                       midpoint = 0, limit = c(0,0.4), space = "Lab", 
#                       name="Fst") +
#  theme_minimal()+ 
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                   size = 12, hjust = 1))+
#  coord_fixed()+
#  xlab("")+
#  ylab("")


##############
## Figure 5  #
##############

# TREE

# Get number of populations
#numpop <- length(unique(geography$location))

#Define OTUs
specieslist <- split(genogeo$sample.id, genogeo$Species2)
#country_list <- split(geography$name, geography$country)
#location_list <- split(geography$name, geography$location)
#ghff_samplesids <- genogeo[genogeo$Species2=="GHFF",]$sample.id
sp_raxml <- read.tree(species_tree)
# Root the tree
#sp_raxml <- root(sp_raxml, which(sp_raxml$tip.label %in% lrff_samplesids))

#Assign OTU
mytree <- groupOTU(sp_raxml, specieslist)

ggtree(mytree, aes(color=group),branch.length="none", layout="rectangular") + 
  scale_color_manual(values=c(rainbow_hcl(6))) + 
  theme(legend.position="right")+ 
  geom_tiplab(size =0.5)+
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1))  + 
  ggplot2::xlim(0, 60)
ggsave("allpops_ml_tree_phylogram_colorbyspecies.pdf")

ggtree(mytree, aes(color=group), layout="rectangular") + geom_tiplab(size=0.5) +
  scale_color_manual(values=c(rainbow_hcl(6))) + 
  theme(legend.position="right")+ 
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label))& as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1)) 
ggsave("allpops_ml_tree_branchlen_colorbyspecies.pdf")


bffsff_list <- split(genogeo$sample.id, genogeo$Broad_region)
bffsff_list <- split(bffsff_popfile$V1, bffsff_popfile$V2)
#country_list <- split(geography$name, geography$country)
#location_list <- split(geography$name, geography$location)
#ghff_samplesids <- genogeo[genogeo$Species2=="GHFF",]$sample.id
bffsff_raxml <- read.tree(bffsff_tree)
# Root the tree
#sp_raxml <- root(sp_raxml, which(sp_raxml$tip.label %in% lrff_samplesids))

#Assign OTU
bffsff_tree <- groupOTU(bffsff_raxml, bffsff_list)
bffsff_cols <- c("grey","#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17")
ggtree(bffsff_tree, aes(color=group),branch.length="none", layout="rectangular") + 
  scale_color_manual(values=c(bffsff_cols )) + 
  theme(legend.position="right")+ 
  geom_tiplab(size =0.75)+
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1))  + 
  ggplot2::xlim(0, 60)
ggsave("allpops_ml_tree_phylogram_colorbyspecies.pdf")

ggtree(bffsff_tree, aes(color=group), layout="rectangular") + geom_tiplab(size=0.5) +
  scale_color_manual(values=c(rainbow_hcl(8))) + 
  theme(legend.position="right")+ 
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label))& as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1)) 
#ggsave("allpops_ml_tree_branchlen_colorbyspecies.pdf")


##############
## Figure 6  #
##############
library(geosphere)
library(adegenet)
library(vegan)
library(ape)
library(vcfR)
library(rdiversity)
library(dplyr)
library(tidyr)
library(ggplot2)
# Ring species figure
genogeo_bffsff <- genogeo_bffsff %>% rowwise() %>% mutate(km_distance_from_sawu=distm(c(121.9167,-10.4833), c(Lon,Lat), fun = distHaversine)/1000)

bffsff_geodist <- genogeo_bffsff %>% select(sample.id,Lon,Lat,Species2,Roost,V2)
bffsff_kinship <- read.delim("2_SNP_Stats/BFFSFF_custom/BFFSFF_plink2_king_kinship.tsv",sep='\t')
bffsff_kinship <- merge(bffsff_kinship,bffsff_geodist, by.x = "X.IID1",by.y = "sample.id")
bffsff_kinship <- merge(bffsff_kinship,bffsff_geodist, by.x = "IID2",by.y = "sample.id")
bffsff_kinship <- bffsff_kinship %>% rowwise() %>% mutate(km_distance=distm(c(Lon.x,Lat.x), c(Lon.y,Lat.y), fun = distHaversine)/1000)

  
bffsff_kinship_exibff <- bffsff_kinship[bffsff_kinship$Species2.x!="IBFF" & bffsff_kinship$Species2.y!="IBFF",]
bffsff_kinship_bff <- bffsff_kinship[bffsff_kinship$Species2.x=="BFF" & bffsff_kinship$Species2.y=="BFF",]
bffsff_kinship_bff <- bffsff_kinship_bff %>% mutate(indo=ifelse((V2.x == "INDO" & V2.y !="INDO")| (V2.y == "INDO" & V2.x !="INDO"),"Indonesia-Australia","Australia-Australia"))

#bffsff_kinship_bff_xpop <- bffsff_kinship[bffsff_kinship$V2.x!=bffsff_kinship$V2.y,]

ggplot(genogeo_bffsff[genogeo_bffsff$Species2!="IBFF",],aes(km_distance_from_sawu,EV1,fill=Species2,shape=Species2)) +
  geom_point(size=4,alpha=0.5)+
  scale_fill_manual(values=c("#E495A5","#C29DDE"))+
  scale_shape_manual(values=c(21,22))+
  ylab("PC1")+
  xlab("Geographic distance from Sawu (km)")+
  theme_bw()+
  theme(text = element_text(size=18),legend.title = element_blank())
ggsave("Distance_from_Sawu_by_PC1.pdf")

ggplot(bffsff_kinship_bff,aes(km_distance,KINSHIP,color=indo,fill=indo)) + geom_point(shape=21,size=2,alpha=0.5)+
  geom_smooth(method='lm', formula= y~x,color="black")+
  scale_fill_manual(values=c("darkred","darkblue"))+
  scale_color_manual(values=c("darkred","darkblue"))+
  ylab("Kinship coefficient")+
  xlab("Geographic distance (km)")+
  theme_bw()+
  theme(text = element_text(size=18),legend.title = element_blank())
ggsave("Kinship_by_distance.pdf")

cor(bffsff_kinship_bff$KINSHIP,bffsff_kinship_bff$km_distance)

cor.test(bffsff_kinship_bff[bffsff_kinship_bff$indo=="Australia-Australia",]$km_distance,
    bffsff_kinship_bff[bffsff_kinship_bff$indo=="Australia-Australia",]$KINSHIP)

cor.test(bffsff_kinship_bff[bffsff_kinship_bff$indo!="Australia-Australia",]$km_distance,
    bffsff_kinship_bff[bffsff_kinship_bff$indo!="Australia-Australia",]$KINSHIP)


bffsff_kinship_bff[bffsff_kinship_bff$KINSHIP> -0.2 & bffsff_kinship_bff$km_distance>1000,] %>% 
  select(V2.x,V2.y) %>% group_by_all %>% count

bffsff_kinship_bff[bffsff_kinship_bff$KINSHIP< -0.2 & bffsff_kinship_bff$km_distance>1000,] %>% 
  select(V2.x,V2.y) %>% group_by_all %>% count

## Mantel test 

bff_km_mat <- bffsff_kinship_bff %>% select(X.IID1,IID2,km_distance)
samples <- unique(c(bff_km_mat$X.IID1, bff_km_mat$IID2))
distance_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples))
rownames(distance_matrix) <- colnames(distance_matrix) <- samples
for (i in 1:nrow(bff_km_mat)) {
  sample1 <- bff_km_mat$X.IID1[i]
  sample2 <- bff_km_mat$IID2[i]
  value <- bff_km_mat$km_distance[i]
  
  # Assign the value to both positions (sample1, sample2) and (sample2, sample1)
  distance_matrix[sample1, sample2] <- value
  distance_matrix[sample2, sample1] <- value
}

bff_kin_mat <- bffsff_kinship_bff %>% select(X.IID1,IID2,KINSHIP)
samples <- unique(c(bff_kin_mat$X.IID1, bff_kin_mat$IID2))
genetic_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples))
rownames(genetic_matrix) <- colnames(genetic_matrix) <- samples
for (i in 1:nrow(bff_kin_mat)) {
  sample1 <- bff_kin_mat$X.IID1[i]
  sample2 <- bff_kin_mat$IID2[i]
  value <- bff_kin_mat$KINSHIP[i]
  
  # Assign the value to both positions (sample1, sample2) and (sample2, sample1)
  genetic_matrix[sample1, sample2] <- value
  genetic_matrix[sample2, sample1] <- value
}

mantel <- mantel(distance_matrix,genetic_matrix)
mantel


# REDO WITHOUT AUSTRALIA-AUSTRALIA
indo_samples <- genogeo_bffsff[genogeo_bffsff$V2=="INDO",]$sample.id

distance_matrix_indo <- matrix(NA, nrow = length(samples), ncol = length(samples))
rownames(distance_matrix_indo) <- colnames(distance_matrix_indo) <- samples
for (i in 1:nrow(bff_km_mat)) {
  sample1 <- bff_km_mat$X.IID1[i]
  sample2 <- bff_km_mat$IID2[i]
  value <- bff_km_mat$km_distance[i]
  if (sample1 %in% indo_samples | sample2 %in% indo_samples) {
  # Assign the value to both positions (sample1, sample2) and (sample2, sample1)
  distance_matrix_indo[sample1, sample2] <- value
  distance_matrix_indo[sample2, sample1] <- value
  } else {
    distance_matrix_indo[sample1, sample2] <- NA
    distance_matrix_indo[sample2, sample1] <- NA
  }
}


genetic_matrix_indo <- matrix(NA, nrow = length(samples), ncol = length(samples))
rownames(genetic_matrix_indo) <- colnames(genetic_matrix_indo) <- samples
for (i in 1:nrow(bff_kin_mat)) {
  sample1 <- bff_kin_mat$X.IID1[i]
  sample2 <- bff_kin_mat$IID2[i]
  value <- bff_kin_mat$KINSHIP[i]
  
  if (sample1 %in% indo_samples | sample2 %in% indo_samples) {
    # Assign the value to both positions (sample1, sample2) and (sample2, sample1)
    genetic_matrix_indo[sample1, sample2] <- value
    genetic_matrix_indo[sample2, sample1] <- value
  } else {
    genetic_matrix_indo[sample1, sample2] <- NA
    genetic_matrix_indo[sample2, sample1] <- NA
  }
}

mantel_indo <- mantel(distance_matrix_indo,genetic_matrix_indo,na.rm=TRUE)
mantel_indo



################
# MISC FIGURES #
################

# Latitude versus PC

g1 <- ggplot(genogeo[genogeo$country=="Australia",],aes(x=Lon,y=EV1,color=broad_region,group=country,shape=species, label=sample.id))+ 
  geom_point() + 
  theme + 
  scale_color_manual(values=c(rainbow_hcl(3))) + 
  #geom_text_repel() + 
  xlab("Longitude") + 
  ylab(ev1pc)+
  geom_smooth(method='lm')
g2 <- ggplot(genogeo[genogeo$country=="Australia",],aes(x=Lon,y=EV2,color=broad_region,group=country,shape=species, label=sample.id))+ 
  geom_point() + 
  theme + 
  scale_color_manual(values=c(rainbow_hcl(3))) + 
  #geom_text_repel() + 
  xlab("Longitude") + 
  ylab(ev2pc)+
  geom_smooth(method='lm')
g3 <- ggplot(genogeo[genogeo$country=="Australia",],aes(x=Lon,y=EV3,color=broad_region,group=country,shape=species, label=sample.id))+ 
  geom_point() + 
  theme + 
  scale_color_manual(values=c(rainbow_hcl(3))) + 
  #geom_text_repel() + 
  xlab("Longitude") + 
  ylab(ev3pc)+
  geom_smooth(method='lm')
g4 <- ggplot(genogeo[genogeo$country=="Australia",],aes(x=Lon,y=EV4,color=broad_region,group=country,shape=species, label=sample.id))+ 
  geom_point() + 
  theme + 
  scale_color_manual(values=c(rainbow_hcl(3))) + 
  #geom_text_repel() + 
  xlab("Longitude") + 
  ylab(ev4pc)+
  geom_smooth(method='lm')

plot_grid(g1,g2,g3,g4)
ggsave("PCA_by_longitude_Australia.pdf")


# Store table

write.table(genogeo,"bff_sff_geno_geo_no_outliers.tsv",sep = '\t',row.names = F,quote = F)

# Highlight putative hybrid samples

hybrids <- c("1510019_LRFF", # Admixed between two GHFF pops
             "1815076BFF_GHFF", # Admixed between BFFxSFF pops
             "1610014_LRFF") # Admixed BFFxSFF

