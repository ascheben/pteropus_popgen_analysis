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
library(geosphere)
library(adegenet)
library(vegan)
library(vcfR)
library(rdiversity)



ff_vcf <- "../../1_SNP_Calling/data/pteropus.vcf.gz"
bffsff_vcf <- "../../1_SNP_Calling/data/palecto_pconspicillatus_broad.vcf.gz"
bffsff_kinship <- read.delim("../../2_SNP_Stats/data/palecto_pconspicillatus.kinship.tsv",sep='\t')

# Read tab-separated popmap ("sample\tgroup") matching names in VCF
popfile <- read.csv("../../0_Metadata/data/pteropus.popmap", header = FALSE, sep = "\t")
geography <- read.csv("../../0_Metadata/data/pteropus_metadata.txt", header = TRUE, sep = "\t") 
bffsff_popfile <- read.csv("../../0_Metadata/data/palecto_pconspicillatus.popmap", header = FALSE, sep = "\t")

# Store tree locations

#species_tree <- "../data/raxml/RAxML_bipartitions.pteropus"
species_tree <- "../data/raxml/RAxML_rerooted.pteropus"
bffsff_tree <- "../data/raxml/RAxML_bipartitions.palecto_pconspicillatus"

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

genogeo <- merge(tabpops,geography,by.x = c("sample.id"),by.y = c("Sample.Identifier"))
genogeo <- genogeo %>% mutate(Species2 = ifelse(species == "BFF" & Country == "Indonesia" & EV1>0,"IBFF",species))

#options(ggrepel.max.overlaps = 100)
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
#ggsave("PCA_by_species_7by8.pdf",height=7,width=8)
ggsave("../figures/PCA_all_species.pdf",height=7,width=8)

genogeo_tab <- genogeo %>% select(sample.id,EV1,EV2,EV3,EV4)
write.table(genogeo_tab,"../figures/Species_PCA.tsv",quote = F,row.names = F,sep = '\t')

##############
## Figure 2  #
##############

# Plot geography

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

#bffsff_cols <- c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17")

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

##############
## Figure 3  #
##############

# Use only BFF and SFF
# Add popmap header
colnames(bffsff_popfile) <- c("sample", "group")
# Set output file names
outname <- "out_bff_sff"
gdsfile <- paste(outname, "gds", sep=".")
# Convert VCF to GDS
#snpgdsVCF2GDS(bffsff_vcf,gdsfile, method="biallelic.only")
# Read in genotypes
genofile <- snpgdsOpen(gdsfile, allow.duplicate=TRUE)
snpset <- snpgdsLDpruning(genofile,autosome.only=FALSE, slide.max.bp = 1)

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
tabpops_bffsff <- merge(tab,bffsff_popfile, by.x = "sample.id" , by.y = "sample" )

genogeo_bffsff <- merge(tabpops_bffsff,geography,by.x = c("sample.id"), by.y = c("Sample.Identifier") )

#############################################################################################################
### WARNING : EV2 THRESHOLD MUST BE MANUALLY CHECKED USING PLOT BELOW TO ENSURE IT ONLY DETECTS TRUE IBFF ###
#############################################################################################################

genogeo_bffsff <- genogeo_bffsff %>% mutate(Species2 = ifelse(Species == "BFF" & Country == "Indonesia" & EV2< -0.1,"IBFF",Species))
genogeo_bffsff <- merge(genogeo_bffsff,bffsff_popfile, by.x = c("sample.id"),by.y = c("sample"))

  
#bffsff_cols <- c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17")

ggplot(genogeo_bffsff,aes(x=EV1,y=EV3,fill=BFF_SFF_Population_ID,shape=Species2, label=sample.id))+ 
  scale_shape_manual(values = c(21,25,22))+ 
  scale_fill_manual(values=c(c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17"))) + 
  geom_point(size=4,alpha=0.5) + 
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  #geom_text_repel() + 
  theme_bw()+
  theme(text = element_text(size=18)) + 
  xlab(ev1pc) + 
  ylab(ev3pc)
ggsave("../figures/PCA_PC1_PC3_bff_sff_species_regions.pdf",height=7,width=10)

ggplot(genogeo_bffsff,aes(x=EV1,y=EV2,fill=BFF_SFF_Population_ID,shape=Species2, label=sample.id))+ 
  scale_shape_manual(values = c(21,25,22))+ 
  scale_fill_manual(values=c(c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17"))) + 
  geom_point(size=4,alpha=0.5) + 
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  #geom_text_repel() + 
  theme_bw()+
  theme(text = element_text(size=18)) + 
  xlab(ev1pc) + 
  ylab(ev2pc)
ggsave("../figures/PCA_PC1_PC2_bff_sff_species_regions.pdf",height=7,width=10)

genogeo_bffsff_tab <- genogeo_bffsff %>% select(sample.id,EV1,EV2,EV3,EV4)
write.table(genogeo_bffsff_tab,"../figures/PCA_bff_sff.tsv",quote = F,row.names = F,sep = '\t')

max.overlaps <- 10
ggplot(genogeo_bffsff,aes(x=EV1,y=EV3,fill=BFF_SFF_Population_ID,shape=Species2, label=sample.id))+ 
  scale_shape_manual(values = c(21,25,22))+ 
  scale_fill_manual(values=c(c("#65BC8C","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17"))) + 
  geom_point(size=4,alpha=0.5) + 
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  geom_text_repel(size=2) + 
  theme_bw()+
  theme(text = element_text(size=18)) + 
  xlab(ev1pc) + 
  ylab(ev3pc)
ggsave("../figures/PCA_PC1_PC3_bff_sff_species_regions_label.pdf",height=7,width=12)

##############
## Figure 5  #
##############

# TREE

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
ggsave("../figures/Phylogram_all_species.pdf",height=10,width=8)

# LRFF and GHFF tree

loc2regiondf <- data.frame(Roost = c("Canungra","Adelaide","Katherine","Tolga","Isa","Wingham","Pt_Macquarie","Casino","Maclean","Townsville","Regents_Park","Charters","Broome","Lismore","Kyogle","Bellingen","Bundgeam","Ingham","Pt_Douglas"), 
                           custom_region = c("ECOAST","SAUS","N_AUS","Wet_tropics","NQ","ECOAST","ECOAST","ECOAST","ECOAST","NQ","ECOAST","NQ","NAUS","ECOAST","ECOAST","ECOAST","ECOAST","NQ","Wet_tropics"))

genogeo <- merge(genogeo,loc2regiondf,by = "Roost", all.x =TRUE)

specieslist <- split(genogeo$sample.id, genogeo$custom_region)
mytreer <- groupOTU(sp_raxml, specieslist)

p <- ggtree(mytreer, aes(color=group), layout="rectangular",branch.length="none") + geom_tiplab(size=1) +
  scale_color_manual(values=c(rainbow_hcl(7))) + 
  theme(legend.position="right")+ 
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label))& as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1)) 
p + geom_text(aes(label=node), hjust=-.3)

p2 <- p %>% collapse(node=246) %>% collapse(node=340) %>% collapse(node=255)
p2
ggsave("../figures/Phylogram_lrff_ghff.pdf",height=10,width=8)


bffsff_list <- split(genogeo$sample.id, genogeo$BFF_SFF_Population_ID)
bffsff_list <- split(bffsff_popfile$sample, bffsff_popfile$group)
bffsff_raxml <- read.tree(bffsff_tree)

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
ggsave("../figures/Phylogram_bff_sff.pdf",height=10,width=8)



##############
## Figure 6  #
##############

genogeo_bffsff <- genogeo_bffsff %>% rowwise() %>% mutate(km_distance_from_sawu=distm(c(121.9167,-10.4833), c(Longitude,Latitude), fun = distHaversine)/1000)
bffsff_geodist <- genogeo_bffsff %>% select(sample.id,Longitude,Latitude,Species2,Roost,BFF_SFF_Population_ID)
bffsff_kinship <- merge(bffsff_kinship,bffsff_geodist, by.x = "X.IID1",by.y = "sample.id")
bffsff_kinship <- merge(bffsff_kinship,bffsff_geodist, by.x = "IID2",by.y = "sample.id")
bffsff_kinship <- bffsff_kinship %>% rowwise() %>% mutate(km_distance=distm(c(Longitude.x,Latitude.x), c(Longitude.y,Latitude.y), fun = distHaversine)/1000)

  
bffsff_kinship_exibff <- bffsff_kinship[bffsff_kinship$Species2.x!="IBFF" & bffsff_kinship$Species2.y!="IBFF",]
bffsff_kinship_bff <- bffsff_kinship[bffsff_kinship$Species2.x=="BFF" & bffsff_kinship$Species2.y=="BFF",]
bffsff_kinship_bff <- bffsff_kinship_bff %>% mutate(indo=ifelse((BFF_SFF_Population_ID.x == "INDO" & BFF_SFF_Population_ID.y !="INDO")| (BFF_SFF_Population_ID.y == "INDO" & BFF_SFF_Population_ID.x !="INDO"),"Indonesia-Australia","Australia-Australia"))

ggplot(genogeo_bffsff[genogeo_bffsff$Species2!="IBFF",],aes(km_distance_from_sawu,EV1,fill=Species2,shape=Species2)) +
  geom_point(size=4,alpha=0.5)+
  scale_fill_manual(values=c("#E495A5","#C29DDE"))+
  scale_shape_manual(values=c(21,22))+
  ylab("PC1")+
  xlab("Geographic distance from Sawu (km)")+
  theme_bw()+
  theme(text = element_text(size=18),legend.title = element_blank())
ggsave("../figures/Distance_from_Sawu_by_PC1.pdf")

ggplot(bffsff_kinship_bff,aes(km_distance,KINSHIP,color=indo,fill=indo)) + geom_point(shape=21,size=2,alpha=0.5)+
  geom_smooth(method='lm', formula= y~x,color="black")+
  scale_fill_manual(values=c("darkred","darkblue"))+
  scale_color_manual(values=c("darkred","darkblue"))+
  ylab("Kinship coefficient")+
  xlab("Geographic distance (km)")+
  theme_bw()+
  theme(text = element_text(size=18),legend.title = element_blank())
ggsave("../figures/Kinship_by_distance.pdf")

#cor(bffsff_kinship_bff$KINSHIP,bffsff_kinship_bff$km_distance)
#cor.test(bffsff_kinship_bff[bffsff_kinship_bff$indo=="Australia-Australia",]$km_distance,
#    bffsff_kinship_bff[bffsff_kinship_bff$indo=="Australia-Australia",]$KINSHIP)
#cor.test(bffsff_kinship_bff[bffsff_kinship_bff$indo!="Australia-Australia",]$km_distance,
#    bffsff_kinship_bff[bffsff_kinship_bff$indo!="Australia-Australia",]$KINSHIP)

#bffsff_kinship_bff[bffsff_kinship_bff$KINSHIP> -0.2 & bffsff_kinship_bff$km_distance>1000,] %>% 
#  select(BFF_SFF_Population_ID.x,BFF_SFF_Population_ID.y) %>% group_by_all %>% count

#bffsff_kinship_bff[bffsff_kinship_bff$KINSHIP< -0.2 & bffsff_kinship_bff$km_distance>1000,] %>% 
#  select(BFF_SFF_Population_ID.x,BFF_SFF_Population_ID.y) %>% group_by_all %>% count

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

mantel_out <- data.frame(
  parameter = c("mantel_pval", "mantel_statistic"),
  value   = c(mantel$signif,mantel$statistic)
)
write.table(mantel_out,"../figures/mantel_test.csv",row.names = F,quote = F,sep=',')

# REDO WITHOUT AUSTRALIA-AUSTRALIA
# indo_samples <- genogeo_bffsff[genogeo_bffsff$BFF_SFF_Population_ID=="INDO",]$sample.id
# 
# distance_matrix_indo <- matrix(NA, nrow = length(samples), ncol = length(samples))
# rownames(distance_matrix_indo) <- colnames(distance_matrix_indo) <- samples
# for (i in 1:nrow(bff_km_mat)) {
#   sample1 <- bff_km_mat$X.IID1[i]
#   sample2 <- bff_km_mat$IID2[i]
#   value <- bff_km_mat$km_distance[i]
#   if (sample1 %in% indo_samples | sample2 %in% indo_samples) {
#   # Assign the value to both positions (sample1, sample2) and (sample2, sample1)
#   distance_matrix_indo[sample1, sample2] <- value
#   distance_matrix_indo[sample2, sample1] <- value
#   } else {
#     distance_matrix_indo[sample1, sample2] <- NA
#     distance_matrix_indo[sample2, sample1] <- NA
#   }
# }
# 
# 
# genetic_matrix_indo <- matrix(NA, nrow = length(samples), ncol = length(samples))
# rownames(genetic_matrix_indo) <- colnames(genetic_matrix_indo) <- samples
# for (i in 1:nrow(bff_kin_mat)) {
#   sample1 <- bff_kin_mat$X.IID1[i]
#   sample2 <- bff_kin_mat$IID2[i]
#   value <- bff_kin_mat$KINSHIP[i]
#   
#   if (sample1 %in% indo_samples | sample2 %in% indo_samples) {
#     # Assign the value to both positions (sample1, sample2) and (sample2, sample1)
#     genetic_matrix_indo[sample1, sample2] <- value
#     genetic_matrix_indo[sample2, sample1] <- value
#   } else {
#     genetic_matrix_indo[sample1, sample2] <- NA
#     genetic_matrix_indo[sample2, sample1] <- NA
#   }
# }
# 
# mantel_indo <- mantel(distance_matrix_indo,genetic_matrix_indo,na.rm=TRUE)
# mantel_indo
