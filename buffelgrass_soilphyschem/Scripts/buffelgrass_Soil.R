# Soil data analysis for buffelgrass project 
# Author: Gabriela Iñigo
# Date: March 2023


# set working directory
setwd("/Users/gabrielainigo/Documents/GitHub/")


# LOAD DATA

# load metadata
metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)

# load soil analysis data 
soil_data = read.csv("buffelgrass_soilphyschem/Data/SoilAnalysisData.csv", header = T)

# add a "Sample" column in metadata for joining it with soil data
metadata$Sample <- row.names(metadata)

# join metadata and soil data
library(dplyr) # for left_join function
soil_metadata <- left_join(soil_data, metadata)

# remove one of CPV2 (row 17)
soil_metadata.fil <- soil_metadata[-c(17),]

# change row names for the values under "Sample" column
row.names(soil_metadata.fil) <- soil_metadata.fil$Sample


# remove one of CPV2 (row 17)
soil_data.fil <- soil_data[-c(17),]

# change row names for the values under "Sample" column
row.names(soil_data.fil) <- soil_data.fil$Sample

# remove "Sample" column
soil_data.sub = subset(soil_data.fil, select = -c(Sample))


# Spearman correlogram
library(ggcorrplot)
soil_corr<-cor(soil_data.sub, method="spearman")
soil_p.mat<-cor_pmat(soil_data.sub)

ggcorrplot(soil_corr, method = "circle", type = "upper", outline.color = "white",
           colors = c("firebrick2", "white", "dodgerblue2"), p.mat = soil_p.mat, insig = "blank")

# Save plot as a pdf
ggsave("Corr_soil.pdf", device = "pdf", width = 12, height = 10, units = "in")


# Eucledean distance matrix
library(vegan)
soil_data.euc <- vegdist(decostand(soil_data.sub, method = "standardize"), 
                        method = "euclidean")

# PERMANOVA 
adonis2(soil_data.euc~SiteInvasion*Cover, data = soil_metadata.fil, permutations = 9999)


# PCA
soil_data.pca <- rda(soil_data.sub, scale=T)
summary(soil_data.pca) #39.5% / 26.6%

soil_pca_metadata<-cbind(soil_metadata.fil, scores(soil_data.pca)$sites)
soil_pca_vectors<-data.frame(pc1 = scores(soil_data.pca)$species[,1], 
                             pc2 = scores(soil_data.pca)$species[,2])

colors_cover<-c("Open_Bare"="burlywood", "Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

library(ggforce)
ggplot(data=soil_pca_metadata, aes(x=PC1, y=PC2, color=Cover))+
  geom_segment(data=soil_pca_vectors, aes(x=0, y=0, xend=pc1*1.5, yend=pc2*1.5), arrow=arrow(length=unit(1/2, "picas")), color="black")+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_hull(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10)+
  annotate("text", x=soil_pca_vectors$pc1*1.6, y=soil_pca_vectors$pc2*1.6, label=rownames(soil_pca_vectors))+
  scale_color_manual(values=colors_cover)+
  labs(x="PC1 (39.5%)", y="PC2 (26.6%)")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("PCA_soil.pdf", device = "pdf", width = 12, height = 10, units = "in")


arid_pca<-rda(arid[,1:18], scale=T)summary(arid_pca) #58.9 / 16.7arid <- cbind(arid, data.frame(pc1 = scores(arid_pca)$sites[,1],                               pc2 = scores(arid_pca)$sites[,2])) arid_pca_species <- data.frame(pc1 = scores(arid_pca)$species[,1],                               pc2 = scores(arid_pca)$species[,2])arid_pcaplot<-ggplot(arid, aes(x = pc1, y = pc2)) +  geom_segment(data = arid_pca_species, aes(x = 0, y = 0, xend = pc1*2, yend = pc2*2), arrow = arrow(length = unit(1/2, "picas")), color = "black", alpha=0.9) +  geom_point(aes(color=Vegetated, fill = Vegetated, size=pH), shape=21, alpha=0.75) +  scale_fill_manual(values=c("darkorange3", "chartreuse3"))+  scale_color_manual(values=c("darkorange3", "chartreuse3"))+  scale_size(range=c(3, 7))+  annotate("text", x = arid_pca_species$pc1*2.2, y = arid_pca_species$pc2*2.2,           label = rownames(arid_pca_species))+  labs(x = "PC1 (58.9%)", y = "PC2 (16.7%)")+  geom_hline(yintercept=0, linetype="dashed")+  geom_vline(xintercept=0, linetype="dashed")+  theme_classic(base_size=15)+  theme(legend.position = "bottom")+  ggtitle("Sonoran desert (metagenomes)")

# PCA with elipses
soil_data.pca <- rda(soil_data.sub, scale=T)
summary(soil_data.pca) #39.5% / 26.6%

soil_pca_metadata<-cbind(soil_metadata.fil, scores(soil_data.pca)$sites)
soil_pca_vectors<-data.frame(pc1 = scores(soil_data.pca)$species[,1], 
                             pc2 = scores(soil_data.pca)$species[,2])

colors_cover<-c("Open_Bare"="burlywood", "Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

library(ggforce)
ggplot(data=soil_pca_metadata, aes(x=PC1, y=PC2, color=Cover))+
  geom_segment(data=soil_pca_vectors, aes(x=0, y=0, xend=pc1*1.5, yend=pc2*1.5), arrow=arrow(length=unit(1/2, "picas")), color="black")+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_polygon(stat = "ellipse", aes(group=SiteInvasion, color="black"), fill=NA) +
  annotate("text", x=soil_pca_vectors$pc1*1.6, y=soil_pca_vectors$pc2*1.6, label=rownames(soil_pca_vectors))+
  scale_color_manual(values=colors_cover)+
  labs(x="PC1 (39.5%)", y="PC2 (26.6%)")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("PCA_ellipses_soil.pdf", device = "pdf", width = 12, height = 10, units = "in")

# PCA with 95% confident interval ellipses
library(ggforce)
ggplot(data=soil_pca_metadata, aes(x=PC1, y=PC2, color=Cover))+
  geom_segment(data=soil_pca_vectors, aes(x=0, y=0, xend=pc1*1.5, yend=pc2*1.5), arrow=arrow(length=unit(1/2, "picas")), color="black")+
  geom_point(aes(shape=SiteInvasion), size=4) +
  stat_ellipse(geom = "polygon", aes(group=SiteInvasion, color="black"), fill=NA) + 
  annotate("text", x=soil_pca_vectors$pc1*1.6, y=soil_pca_vectors$pc2*1.6, label=rownames(soil_pca_vectors))+
  scale_color_manual(values=colors_cover)+
  labs(x="PC1 (39.5%)", y="PC2 (26.6%)")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("PCA_ellipses_95_soil.pdf", device = "pdf", width = 12, height = 10, units = "in")

# Linear models

# linear model with water content as response variable
soil_WaterCon.lm <- lm(WaterContent ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_WaterCon.lm)

# linear model with pH as response variable
soil_pH.lm <- lm(pH ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_pH.lm)

# linear model with EC as response variable
soil_EC.lm <- lm(EC ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_EC.lm)

# linear model with N % as response variable
soil_N.lm <- lm(N ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_N.lm)

# linear model with C % as response variable
soil_C.lm <- lm(C ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_C.lm)

# linear model with Fe as response variable
soil_Fe.lm <- lm(Fe ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_Fe.lm)

# linear model with Cu as response variable
soil_Cu.lm <- lm(Cu ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_Cu.lm)

# linear model with Zn as response variable
soil_Zn.lm <- lm(Zn ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_Zn.lm)

# linear model with K as response variable
soil_K.lm <- lm(K ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_K.lm)

# linear model with Mg as response variable
soil_Mg.lm <- lm(Mg ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_Mg.lm)

# linear model with Ca as response variable
soil_Ca.lm <- lm(Ca ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_Ca.lm)

# linear model with P as response variable
soil_P.lm <- lm(P ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_P.lm)

# linear model with Mn as response variable
soil_Mn.lm <- lm(Mn ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_Mn.lm)

# linear model with S as response variable
soil_S.lm <- lm(S ~ SiteInvasion * Cover, data = soil_metadata.fil)
anova(soil_S.lm)
