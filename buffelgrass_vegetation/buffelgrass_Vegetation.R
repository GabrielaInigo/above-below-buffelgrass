# Vegetation data analysis for buffelgrass project
# Author: Gabriela Iñigo
# Date: December 2023

# set working directory
setwd("/Users/gabrielainigo/Documents/GitHub/above-below-buffelgrass")

# LOAD DATA

# load vegetation data
veg <- read.csv("/Users/gabrielainigo/Documents/GitHub/above-below-buffelgrass/buffelgrass_vegetation/VegAbundanceTable.csv", row.names = 1)

# remove first and second column of veg
veg <- veg[,c(-1,-2)]

# load metadata
metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)

# Create a metadata without opens
#veg_metadata <- metadata[!metadata$Vegetation %in% c("Open Bare", "Open Native", "Open Buffel),] 

# create a metadata without "open bare" since the veg data just has "open native" at the natural site and "open buffel" at the transformed site
veg_metadata1 <- metadata[!(metadata$Vegetation %in% c("Open Bare")),]

# Remove "Open Native" at transformed site from veg_metadata1
veg_metadata2 <- veg_metadata1[!(veg_metadata1$sampleid %in% c("ion1","ion2","ion3","ion4","ion5")),]

# Remove "Open Buffel" at natural site from veg_metadata2
veg_metadata <- veg_metadata2[!(veg_metadata2$sampleid %in% c("cob1","cob2","cob3","cob4","cob5")),]


# DATA EXPLORATION

# make a histogram of the veg table to know the distribution of the data
hist(rowSums(veg))

# get stat summary of the data
summary(rowSums(veg))

# match rownames of veg with rownames of metadata
veg = veg[rownames(veg_metadata),]

# make sure all are TRUE
rownames(veg)==rownames(veg_metadata)


# RICHNESS & DIVERSITY

# richness
library(vegan) # for specnumber function
veg_metadata$Richness.veg <- specnumber(veg, MARGIN = 1)

# diversity 
veg_metadata$Shannon.veg <- diversity(veg, MARGIN = 1)


# convert Site Invasion and CoverGeneral into factors
veg_metadata$SiteInvasion<-as.factor(veg_metadata$SiteInvasion)
veg_metadata$CoverGeneral<-as.factor(veg_metadata$CoverGeneral)
veg_metadata$CoverGeneral <- factor(veg_metadata$CoverGeneral, levels=c("Open", "Mesquite", "IronWood", "PaloVerde"))
veg_metadata$Cover<-as.factor(veg_metadata$Cover)
veg_metadata$Cover <- factor(veg_metadata$Cover, levels=c("Open_withBuffel", "Open_withNatives",
                                                          "Mesquite_noBuffel", "Mesquite_withBuffel",
                                                          "IronWood_noBuffel", "IronWood_withBuffel",
                                                          "PaloVerde_noBuffel", "PaloVerde_withBuffel"))

colors_cover<-c("Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

# vegetation richness plot
library(ggplot2) #for ggplot function
library(ggh4x) #for nested facet
ggplot(veg_metadata, aes(x=Cover, y=Richness.veg, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Vegetation richness") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

# save plot as a pdf
ggsave("Vegetation_Richness.pdf", device = "pdf", width = 12, height = 10, units = "in")

# vegetation richness linear model
veg_richness_lm<-lm(Richness.veg ~ SiteInvasion * Cover, data = veg_metadata)
anova(veg_richness_lm)

# r2
library(performance) # for r2 function
r2(veg_richness_lm)


# vegetation diversity plot
library(ggplot2) #for ggplot function
library(ggh4x) #for nested facet
ggplot(veg_metadata, aes(x=Cover, y=Shannon.veg, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Vegetation diversity") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

# save plot as a pdf
ggsave("Vegetation_Diversity.pdf", device = "pdf", width = 12, height = 10, units = "in")

# vegetation diversity linear model
veg_diversity_lm<-lm(Shannon.veg ~ SiteInvasion * Cover, data = veg_metadata)
anova(veg_diversity_lm)

# r2
library(performance) # for r2 function
r2(veg_diversity_lm)


# COMPOSITION SIMILARITY

library(vegan)
veg_bray <- vegdist(veg, method="bray")
veg_nmds <- metaMDS(veg_bray, k=2, try = 100)
veg_metadata$Axis01 = veg_nmds$points[,1]
veg_metadata$Axis02 = veg_nmds$points[,2]

veg_nmds$stress # 0.1536455

library(ggforce)
ggplot(veg_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_hull(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10)+
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Vegetation_NMDS.pdf", device = "pdf", width = 12, height = 10, units = "in")

# PERMANOVA
adonis2(veg_bray~SiteInvasion*Cover, data = veg_metadata, permutations = 9999)

# Betadisper
bd <- betadisper(veg_bray, veg_metadata$SiteInvasion) #Avergae distance to median of induced site:0.29 and natural:0.43
# test if weather or not there is significant difference in variation between sites
anova(bd) # p-value: 0.007273

# Envfit
veg_envfit <- envfit(veg_nmds, veg, permutations = 999, na.rm = TRUE)

# run veg_envfit
veg_envfit

# plot NMDS and envfit results in base R 
plot(veg_nmds)
plot(veg_envfit)

# Extract coordinates from veg_enfit and multiply them by ordiArrowMul(veg_envfit) to have them in correct proportion
veg_envfit_coord <- as.data.frame(scores(veg_envfit, "vectors")) * ordiArrowMul(veg_envfit)

# Plot NMDS and envfit
library(ggforce)
ggplot(veg_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_hull(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10)+
  scale_color_manual(values=colors_cover)+
  geom_segment(data = veg_envfit_coord, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = veg_envfit_coord, aes(x=NMDS1, y=NMDS2), colour = "grey30", 
           fontface = "bold", label = colnames(veg)) + 
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Vegetation_Envfit.pdf", device = "pdf", width = 12, height = 10, units = "in")

# Envfit coordinates of only significant species
veg_envfit_coord.fil <- veg_envfit_coord[c("abuabu","ariads","boubar",
                                         "cencil","chanic","euperi","eupflo",
                                         "ipohed","lopsch","setmac",
                                         "steala", "tidlan","bare_soil"),]

# plot nmds and envfit of only significant species
# Plot NMDS and envfit
library(ggforce)
ggplot(veg_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_hull(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10)+
  scale_color_manual(values=colors_cover)+
  geom_segment(data = veg_envfit_coord.fil, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = veg_envfit_coord.fil, aes(x=NMDS1, y=NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(veg_envfit_coord.fil)) + 
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Vegetation_Envfit_sig.pdf", device = "pdf", width = 12, height = 10, units = "in")

# Vegetation NMDS with envfit of only sig species and 95% confince interval ellipses
library(ggforce)
ggplot(veg_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  stat_ellipse(geom = "polygon", aes(group=SiteInvasion, color="black"), fill=NA)+
  scale_color_manual(values=colors_cover)+
  geom_segment(data = veg_envfit_coord.fil, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = veg_envfit_coord.fil, aes(x=NMDS1, y=NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(veg_envfit_coord.fil)) + 
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Vegetation_Envfit_sig_95.pdf", device = "pdf", width = 12, height = 10, units = "in")


# Hierarchical clustering

# bray-curtis distance matrix
veg_bray <- vegdist(veg, method="bray")

# average
veg_clust.a <- hclust(veg_bray, method="average")
plot(veg_clust.a)


# PCA

# add Sample column to veg for joining with veg_metadata
veg$Sample <- row.names(veg)

# add Sample column to veg_metadata
veg_metadata$Sample <- row.names(veg_metadata)

# join veg_metadata and veg
library(dplyr) # for left_join function
veg_data_meta <- left_join(veg, veg_metadata)

# change row names for the values under "Sample" column
row.names(veg_pca) <- veg_pca$Sample

# remove "Sample" column
veg = subset(veg, select = -c(Sample))

# remove "sampleid" column
veg = subset(veg, select = -c(sampleid))

veg.pca <- rda(veg, scale=T)
summary(veg.pca) #10.48% / 7.90%

veg_pca_metadata<-cbind(veg_data_meta, scores(veg.pca)$sites)
veg_pca_vectors<-data.frame(pc1 = scores(veg.pca)$species[,1], 
                             pc2 = scores(veg.pca)$species[,2])

colors_cover<-c("Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

library(ggforce)
ggplot(data=veg_pca_metadata, aes(x=PC1, y=PC2, color=Cover))+
  geom_segment(data=veg_pca_vectors, aes(x=0, y=0, xend=pc1*1.5, yend=pc2*1.5), arrow=arrow(length=unit(1/2, "picas")), color="black")+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_hull(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10)+
  annotate("text", x=veg_pca_vectors$pc1*1.6, y=veg_pca_vectors$pc2*1.6, label=rownames(veg_pca_vectors))+
  scale_color_manual(values=colors_cover)+
  labs(x="PC1 (10.48%)", y="PC2 (7.90%)")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("PCA_veg.pdf", device = "pdf", width = 12, height = 10, units = "in")



# SPECIES ACCUMULATION CURVE

# download needed packages
install.packages("BiodiversityR")
library(BiodiversityR)

install.packages("ggsci")
library(ggsci)

install.packages("readxl")
library(readxl)

# create copy of veg 
veg_2 <- veg

# convert cover percent to presence absence
veg_2[veg_2 > 0] <- 1

# Get the species accumulation result
Accum.1 <- accumcomp(veg_2, y = veg_metadata, factor = 'Cover',
                     method = 'exact', conditioned = FALSE, plotit = FALSE)
Accum.1



# Render species accumulation data in the long format
accum.long1 <- accumcomp.long(Accum.1, ci=NA, label.freq=8)
head(accum.long1)

# Plot the species accumulation curve of both sites
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

library(ggplot2)
ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), linewidth=2) +
  geom_point(data=subset(accum.long1, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE) + 
  BioR.theme +
  scale_colour_npg() +
  labs(x = "Samples", y = "Plant Species", colour = "Cover", shape = "Cover")

# Save plot as a pdf
ggsave("SAC_veg_Cover.pdf", device = "pdf", width = 12, height = 10, units = "in")

# Species accumulation curve of transformed site

# data
veg_2_t <- veg_2[!(row.names(veg_2) %in% c("NPV1","NPV2","NPV3","NPV4","NPV5",
                          "NPF1","NPF2","NPF3","NPF4","NPF5",
                          "NM1","NM2","NM3","NM4","NM5",
                          "NON1","NON2","NON3","NON4","NON5",
                          "CPV1","CPV2","CPV3","CPV4","CPV5",
                          "CPF1","CPF2","CPF3","CPF4","CPF5",
                          "CM1","CM2","CM3","CM4","CM5")),]
# metadata
veg_metadata_t <- veg_metadata[!(veg_metadata$SiteInvasion %in% c("Natural")),]

# Get the species accumulation result
Accum.1.t <- accumcomp(veg_2_t, y = veg_metadata_t, factor = 'Cover',
                     method = 'exact', conditioned = FALSE, plotit = FALSE)

# Render species accumulation data in the long format
accum.long1.t <- accumcomp.long(Accum.1.t, ci=NA, label.freq=4)
head(accum.long1.t)

# Plot the species accumulation curve of the transformed site
library(ggplot2)
ggplot(data=accum.long1.t, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), linewidth=2) +
  geom_point(data=subset(accum.long1.t, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE) + 
  BioR.theme +
  scale_colour_npg() +
  labs(x = "Samples", y = "Plant Species", colour = "Cover", shape = "Cover")

# Save plot as a pdf
ggsave("SAC_veg_Transformed.pdf", device = "pdf", width = 12, height = 10, units = "in")


# Species accumulation curve of Natural site

# data
veg_2_n <- veg_2[!(row.names(veg_2) %in% c("IPV1","IPV2","IPV3","IPV4","IPV5",
                                                "IPF1","IPF2","IPF3","IPF4","IPF5",
                                                "IM1", "IM2","IM3","IM4","IM5",
                                                "IOB1","IOB2","IOB3","IOB4","IOB5")),]

# metadata
veg_metadata_n <- veg_metadata[!(veg_metadata$SiteInvasion %in% c("Induced")),]

# Get the species accumulation result
Accum.1.n <- accumcomp(veg_2_n, y = veg_metadata_n, factor = 'Cover',
                       method = 'exact', conditioned = FALSE, plotit = FALSE)

# Render species accumulation data in the long format
accum.long1.n <- accumcomp.long(Accum.1.n, ci=NA, label.freq=7)
head(accum.long1.n)

# Plot the species accumulation curve of the natural site
library(ggplot2)
ggplot(data=accum.long1.n, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), linewidth=2) +
  geom_point(data=subset(accum.long1.n, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, alpha=0.3))), 
              show.legend=FALSE) + 
  BioR.theme +
  scale_colour_npg() +
  labs(x = "Samples", y = "Plant Species", colour = "Cover", shape = "Cover")

# Save plot as a pdf
ggsave("SAC_veg_Natural.pdf", device = "pdf", width = 12, height = 10, units = "in")


# Get the species accumulation result
Accum.2 <- accumcomp(veg_2, y = veg_metadata, factor = 'SiteInvasion',
                     method = 'exact', conditioned = FALSE, plotit = FALSE)
Accum.2

# Render species accumulation data in the long format
accum.long2 <- accumcomp.long(Accum.2, ci=NA, label.freq=2)
head(accum.long2)

library(ggplot2)
ggplot(data=accum.long2, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), linewidth=2) +
  geom_point(data=subset(accum.long2, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, alpha=0.3))), 
              show.legend=FALSE) + 
  BioR.theme +
  scale_colour_npg() +
  labs(x = "Samples", y = "Plant Species", colour = "SiteInvasion", shape = "SiteInvasion")

# Save plot as a pdf
ggsave("SAC_veg_Sites.pdf", device = "pdf", width = 12, height = 10, units = "in")



# get the proportion of invasive and ruderal species

# create copy of veg 
veg_2 <- veg

# convert cover percent to presence absence
veg_2[veg_2 > 0] <- 1

# sum rows of veg_2 and add the results to veg_metadata into a column called "Totalspecies"
veg_metadata$TotalSpecies <- rowSums(veg_2)

# make copy of veg_2 with only invasive and ruderal species
veg_inv_rud <- veg_2[,c("ambcon","ariads","chlvir","galapa","ipohed","phaspi", "cencil")]

# sum rows of veg_inv_rud and add results to veg_metadat into a column called "RuderalandInvasive"
veg_metadata$RuderalandInvasive <- rowSums(veg_inv_rud)

# calculate the proportion of invasive species per sample and add the result to veg_metadata into a column called "Proportion"
veg_metadata$Proportion <- (veg_metadata$RuderalandInvasive/veg_metadata$TotalSpecies)*100

# proportion of invasive and ruderal species
library(ggplot2) #for ggplot function
library(ggh4x) #for nested facet
ggplot(veg_metadata, aes(x=Cover, y=Proportion, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Proportion of Invasive and Ruderal Species") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

# save plot as a pdf
ggsave("Vegetation_proportion.pdf", device = "pdf", width = 12, height = 10, units = "in")


# proportion of invasive and ruderal species linear model
proportion_lm<-lm(Proportion ~ SiteInvasion * Cover, data = veg_metadata) 
anova(proportion_lm)



# r2
library(performance) # for r2 function
r2(proportion_lm)


# Tree metrics stats

# Load tree metrics data
tree_metrics <- read.csv("/Users/gabrielainigo/Documents/GitHub/buffelgrass_vegetation/TreeMetrics.csv", row.names = 1)

# load metadata
metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)

# Create a metadata without opens
tree_metadata <- metadata[(!metadata$Vegetation %in% c("Open Bare", "Open Native", "Open Buffel")),]

# match rownames of veg with rownames of metadata
tree_metrics = tree_metrics[rownames(tree_metadata),]

# make sure all are TRUE
rownames(tree_metrics)==rownames(tree_metadata)

# Add Height column to tree_metadata
tree_metadata$Height = tree_metrics$Height

# Add Diameter column to tree_metadata
tree_metadata$Diameter = tree_metrics$Diameter

# Add TreeCover column to tree_metadata
tree_metadata$TreeCover = tree_metrics$TreeCover

# Add Ellipsoid column to tree_metadata
tree_metadata$Ellipsoid = tree_metrics$Ellipsoid

# Effect of site, cover and the interaction between them on the height of trees
height_lm<-lm(Height ~ SiteInvasion * Cover, data = tree_metadata)
anova(height_lm) # Not significant

# Effect of site, cover and the interaction between them on the diameter of trees
diameter_lm<-lm(Diameter ~ SiteInvasion * Cover, data = tree_metadata)
anova(diameter_lm) # Not significant

# Effect of site, cover and the interaction between them on the cover of trees
TreeCover_lm<-lm(TreeCover ~ SiteInvasion * Cover, data = tree_metadata)
anova(TreeCover_lm) # Not significant

# Effect of site, cover and the interaction between them on the ellipsoid of trees
Ellipsoid_lm<-lm(Ellipsoid ~ SiteInvasion * Cover, data = tree_metadata)
anova(Ellipsoid_lm) # Not significant
