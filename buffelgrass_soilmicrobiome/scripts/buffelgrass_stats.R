# Data analysis for buffelgrass project 
# Author: Gabriela Iñigo
# Date: March 2023

library(vegan) # multivar stats in ecology
library(data.table)
library(dbplyr) # data wrangling
library(tidyverse)
library(ggplot2) # pretty plots
library(lme4) # mixed effect models
library(lmerTest) # model testing for lme4
library(performance)
library(see)
library(ggforce)
library(modelbased) #estimating means of linear models
library(ggforce) #NMDS polygons

setwd("/Users/gabrielainigo/Documents/GitHub/")

#----------------------------------------------------------------------- 
#Load the data
#----------------------------------------------------------------------- 

metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)

# 16S 
b_taxa = read.csv("buffelgrass_soilmicrobiome/data/16S/taxa_table_16S.csv", header = T)

b_asv = read.csv("buffelgrass_soilmicrobiome/data/16S/asv_table_16S.csv", header = T)

# ITS
f_taxa = read.csv("buffelgrass_soilmicrobiome/data/ITS/taxa_table_ITS.csv", header = T)

f_asv = read.csv("buffelgrass_soilmicrobiome/data/ITS/asv_table_ITS.csv", header = T)

# soil analysis data 
soil_data = read.csv("buffelgrass_soilphyschem/Data/SoilAnalysisData.csv", header = T)

# Check out blanks 16S
b_blanks = b_asv[c(1,2,3,4,5,6,7),]

# Check if there is no contamination by adding the blanks
rowSums(b_blanks)

# remove blanks from asv table
b_asv = data.frame(b_asv[-c(1,2,3,4,5,6,7),])


# Check out blanks ITS
f_blanks = f_asv[c(1,2,3,4,5,6,7),] 

# Check if there is no contamination by adding the blanks
rowSums(f_blanks)

# remove blanks from asv table
f_asv = data.frame(f_asv[-c(1,2,3,4,5,6,7),])

# remove unwanted taxonomy from bacteria
fil_b_taxa = b_taxa[!b_taxa$Order %in% "Chloroplast",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!fil_b_taxa$Family %in% "Mitochondria",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!fil_b_taxa$Kingdom %in% "Eukaryota",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!(is.na(fil_b_taxa$Kingdom)),] %>% droplevels()

# remove from ASV table
fil_b_asv <- b_asv[,rownames(fil_b_taxa)]

# remove unwanted taxonomy from fungi
fil_f_taxa = f_taxa[f_taxa$Kingdom %in% "k__Fungi",] %>% droplevels() 

# remove from ASV table
fil_f_asv <- f_asv[,rownames(fil_f_taxa)]

# change row names of filtered asv tables to upper case
# this is necessary to match rownames of asv tables with rownames of metadata
rownames(fil_b_asv) <- toupper(rownames(fil_b_asv))
rownames(fil_f_asv) <- toupper(rownames(fil_f_asv))


#----------------------------------------------------------------------- 
# Data exploration
#----------------------------------------------------------------------- 

#16S
hist(rowSums(fil_b_asv))
summary(rowSums(fil_b_asv))
sd(rowSums(fil_b_asv))

sort(rowSums(fil_b_asv)) #CM2 is the one with ~4000

# Remove cm2
fil_b_asv  <- fil_b_asv[!(row.names(fil_b_asv) == "CM2"),]
b_metadata <- metadata
b_metadata = b_metadata[rownames(fil_b_asv),] #matches rownames metadata with rownames asv table
rownames(b_metadata)==rownames(fil_b_asv) #make sure all are TRUE

summary(rowSums(fil_b_asv)) # min = 40286

#ITS
hist(rowSums(fil_f_asv)) 
summary(rowSums(fil_f_asv))
sd(rowSums(fil_f_asv))
sort(rowSums(fil_f_asv)) # IOB5 is the one with ~1300

# Remove iob5
fil_f_asv  <- fil_f_asv[!(row.names(fil_f_asv) == "IOB5"),]
f_metadata <- metadata
f_metadata = f_metadata[rownames(fil_f_asv),] #matches rownames metadata with rownames asv table
rownames(f_metadata)==rownames(fil_f_asv) #make sure all are TRUE

summary(rowSums(fil_f_asv))# min 27640

# create copy abundance matrix
fil_asv_ab <- fil_b_asv
# convert abundance to presence absence
fil_asv_ab[fil_asv_ab > 0] <- 1
#calculate the mean of the total number of phylotypes in samples
mean(rowSums(fil_asv_ab))
range(rowSums(fil_asv_ab))
summary(rowSums(fil_asv_ab))

sd(rowSums(fil_asv_ab)) # cosas para reportar en el paper

#ITS
# create copy abundance matrix
fil_asv_ab <- fil_f_asv
# convert abundance to presence abscence
fil_asv_ab[fil_asv_ab > 0] <- 1
#calculate the mean of the total number of phylotypes in samples
mean(rowSums(fil_asv_ab))
range(rowSums(fil_asv_ab))
summary(rowSums(fil_asv_ab))
sd(rowSums(fil_asv_ab))

# Soil 

# add a "Sample" column in metadata for joining it with soil data
metadata$Sample <- row.names(metadata)

# join metadata and soil data
soil_metadata <- left_join(soil_data, metadata)

# remove one of CPV2 (row 17)
soil_metadata.fil <- soil_metadata[-c(17),]

# change row names for the values under "Sample" column
row.names(soil_metadata.fil) <- soil_metadata.fil$Sample

# ---------------------------------------------------------------------------------
## ALPHA DIVERSITY ##
# ---------------------------------------------------------------------------------

# 16S
b_metadata$Richness.16S.rar <- specnumber(rrarefy(fil_b_asv, sample=40000)) # remove small outliers
b_metadata$Shannon.16S <- diversity(rrarefy(fil_b_asv, sample=40000), index = "shannon")


# ITS
f_metadata$Richness.ITS.rar<-specnumber(rrarefy(fil_f_asv, sample=20000))
f_metadata$Shannon.ITS <- diversity(rrarefy(fil_f_asv, sample=20000), index = "shannon")


# ---------------------------------------------------------------------------------
## Soil cover/Vegetation Analysis ##
# ---------------------------------------------------------------------------------

# susbset of the bacteria metadata without trees
cover_b_metadata <- subset.data.frame(b_metadata[c("IOB1","IOB2", "IOB3", "IOB4", "IOB5",
                                           "ION1", "ION2", "ION3", "ION4", "ION5",
                                           "IOD1", "IOD2", "IOD3", "IOD4", "IOD5",
                                           "NON1", "NON2", "NON3", "NON4", "NON5",
                                           "NOD1", "NOD2", "NOD3", "NOD4", "NOD5", 
                                           "COB1", "COB2", "COB3", "COB4", "COB5"),])

# replace site "C" with "N" in order to have 2 sites instead of 3
cover_b_metadata$Site[cover_b_metadata$Site == "C"] <- "N" 

# duplicate the vegetation column and call it "SoilCover"
cover_b_metadata$SoilCover = cover_b_metadata$Vegetation

# susbset of the fungi metadata without trees
cover_f_metadata <- subset.data.frame(f_metadata[c("IOB1","IOB2", "IOB3", "IOB4",
                                                  "ION1", "ION2", "ION3", "ION4", "ION5",
                                                  "IOD1", "IOD2", "IOD3", "IOD4", "IOD5",
                                                  "NON1", "NON2", "NON3", "NON4", "NON5",
                                                  "NOD1", "NOD2", "NOD3", "NOD4", "NOD5", 
                                                  "COB1", "COB2", "COB3", "COB4", "COB5"),])

# replace site "C" with "N" in order to have 2 sites instead of 3
cover_f_metadata$Site[cover_f_metadata$Site == "C"] <- "N" 

# duplicate the vegetation column and call it "SoilCover"
cover_f_metadata$SoilCover = cover_f_metadata$Vegetation

# linear models & plots 

# 16S

# plot of Bacterial and Archaeal richness (observed means and standard erros)
ggplot(data = cover_b_metadata, aes(x = Site, y = Richness.16S.rar, color = SoilCover)) +
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data = "mean_se") +
  geom_jitter(aes(x=Site, y=Richness.16S.rar), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Site") +
  ylab("Bacterial and Archaeal Richness") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB")) +
  scale_x_discrete(labels=c("I" = "Induced", "N" = "Natural Invasion")) +
  facet_grid(. ~ Site, scales = "free", space = "free") +
  theme(strip.text.x = element_blank()) 
#expand_limits(y = 0) this makes the y-axis to start at 0
 

# Save plot as a pdf
ggsave("Rich_bact_cover.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# linear model with bacterial richness as response variable
cover_b_richness.lm <- lm(Richness.16S.rar ~ Site * SoilCover, data = cover_b_metadata)
anova(cover_b_richness.lm)

# estimate means of the previous model to plot them
cover_b_richness.means <- estimate_means(cover_b_richness.lm)
cover_b_richness.means

# plot of Bacterial and Archaeal richness (estimated means and confidence intervals)
ggplot(data = cover_b_richness.means, aes(x = Site, y = Mean, color = SoilCover)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge2(width=1)) +
  geom_jitter(aes(x=Site, y=Mean), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Site") +
  ylab("Bacterial and Archaeal Richness") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB")) +
  scale_x_discrete(labels=c("I" = "Induced", "N" = "Natural Invasion")) +
  facet_grid(. ~ Site, scales = "free", space = "free") +
  theme(strip.text.x = element_blank())

# Save plot as a pdf
ggsave("Rich_bact_cover_lm.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# plot of bacterial and archaeal diversity (observed mean and standard errors)
ggplot(data = cover_b_metadata, aes(x = Site, y = Shannon.16S, color = SoilCover)) +
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data = "mean_se") +
  geom_jitter(aes(x=Site, y=Shannon.16S), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Site") +
  ylab("Bacterial and Archaeal Diversity") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB")) +
  scale_x_discrete(labels=c("I" = "Induced", "N" = "Natural Invasion")) +
  facet_grid(. ~ Site, scales = "free", space = "free") +
  theme(strip.text.x = element_blank())

# save plot as a pdf
ggsave("Div_bact_cover.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# linear model with bacterial diversity as response variable
cover_b_diversity.lm <- lm(Shannon.16S ~ Site * SoilCover, data = cover_b_metadata)
summary(cover_b_diversity.lm)
anova(cover_b_diversity.lm)

# estimate means of the previous model to plot them
cover_b_diversity.means <- estimate_means(cover_b_diversity.lm)
cover_b_diversity.means

# plot of Bacterial and Archaeal richness (estimated means and confidence intervals)
ggplot(data = cover_b_diversity.means, aes(x = Site, y = Mean, color = SoilCover)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge2(width=1)) +
  geom_jitter(aes(x=Site, y=Mean), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Site") +
  ylab("Bacterial and Archaeal Diversity") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB")) +
  scale_x_discrete(labels=c("I" = "Induced", "N" = "Natural Invasion")) +
  facet_grid(. ~ Site, scales = "free", space = "free") +
  theme(strip.text.x = element_blank())

# save plot as a pdf
ggsave("Div_bact_cover_lm.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# ITS

# plot of fungal richness (observed means and standard errors)
ggplot(data = cover_f_metadata, aes(x = Site, y = Richness.ITS.rar, color = SoilCover)) +
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data="mean_se") +
  geom_jitter(aes(x=Site, y=Richness.ITS.rar), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Site") +
  ylab("Fungal Richness") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB")) +
  scale_x_discrete(labels=c("I" = "Induced", "N" = "Natural Invasion")) +
  facet_grid(. ~ Site, scales = "free", space = "free") +
  theme(strip.text.x = element_blank())

# save plot as pdf
ggsave("Rich_fung_cover.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# linear model with fungal richness as response variable
cover_f_richness.lm <- lm(Richness.ITS.rar ~ Site * SoilCover, data = cover_f_metadata)
summary(cover_f_richness.lm)
anova(cover_f_richness.lm)

# estimate means of the previous model to plot them
cover_f_richness.means <- estimate_means(cover_f_richness.lm)
cover_f_richness.means

# plot of fungal richness (estimated means and confidence intervals)
ggplot(data = cover_f_richness.means, aes(x = Site, y = Mean, color = SoilCover)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge2(width=1)) +
  geom_jitter(aes(x=Site, y=Mean), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Site") +
  ylab("Fungal Richness") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB")) +
  scale_x_discrete(labels=c("I" = "Induced", "N" = "Natural Invasion")) +
  facet_grid(. ~ Site, scales = "free", space = "free") +
  theme(strip.text.x = element_blank())

# save plot as pdf
ggsave("Rich_fung_cover_lm.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# plot of fungal diversity (observed means and standard errors)
ggplot(data = cover_f_metadata, aes(x = Site, y = Shannon.ITS, color = SoilCover))+
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data="mean_se") +
  geom_jitter(aes(x=Site, y=Shannon.ITS), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Site") +
  ylab("Fungal Diversity") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB"))+
  scale_x_discrete(labels=c("I" = "Induced", "N" = "Natural Invasion")) +
  facet_grid(. ~ Site, scales = "free", space = "free") +
  theme(strip.text.x = element_blank())

# save plot as pdf
ggsave("Div_fung_cover.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# linear model with fungal diversity as response variable
cover_f_diversity.lm <- lm(Shannon.ITS ~ Site * SoilCover, data = cover_f_metadata)
summary(cover_f_diversity.lm)
anova(cover_f_diversity.lm)

# estimate means of the previous model to plot them
cover_f_diversity.means <- estimate_means(cover_f_diversity.lm)
cover_f_diversity.means

# plot of fungal richness (estimated means and confidence intervals)
ggplot(data = cover_f_diversity.means, aes(x = Site, y = Mean, color = SoilCover)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge2(width=1)) +
  geom_jitter(aes(x=Site, y=Mean), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Site") +
  ylab("Fungal Diversity") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB")) +
  scale_x_discrete(labels=c("I" = "Induced", "N" = "Natural Invasion")) +
  facet_grid(. ~ Site, scales = "free", space = "free") +
  theme(strip.text.x = element_blank())

# save plot as pdf
ggsave("Div_fung_cover_lm.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# Community composition 

# 16S

cover_fil_b_asv <- subset.data.frame(fil_b_asv[c("IOB1","IOB2", "IOB3", "IOB4", "IOB5",
                                                  "ION1", "ION2", "ION3", "ION4", "ION5",
                                                  "IOD1", "IOD2", "IOD3", "IOD4", "IOD5",
                                                  "NON1", "NON2", "NON3", "NON4", "NON5",
                                                  "NOD1", "NOD2", "NOD3", "NOD4", "NOD5", 
                                                  "COB1", "COB2", "COB3", "COB4", "COB5"),])


rownames(cover_fil_b_asv)
cover_b_metadata$sampleid

cover_b_asv_bray <- vegdist(cover_fil_b_asv, method="bray")
cover_b_asv_nmds <- metaMDS(cover_b_asv_bray, k=2, try = 100)
cover_b_metadata$Axis01 = cover_b_asv_nmds$points[,1]
cover_b_metadata$Axis02 = cover_b_asv_nmds$points[,2]

cover_b_asv_nmds$stress

ggplot(cover_b_metadata,aes(Axis01, Axis02))+
  geom_point(aes(color=Site, shape=SoilCover), size=4)+
  labs(title="16S NMDS ordination")+
  theme_classic(base_size = 14)

# perform permanova
adonis2(cover_b_asv_bray ~ SoilCover*Site, data = cover_b_metadata)

# assign labels to use in facet_wrap
# labels_cover <- c("I" = "Induced", "N" = "Natural Invasion")

# plot bacterial NMDS
ggplot(data = cover_b_metadata, aes(Axis01, Axis02)) +
  geom_point(aes(color=SoilCover), size=4) +
  geom_mark_hull(aes(group = Site, label = Site), concavity = 10) +
  labs(title="Bacterial NMDS ordination (Stress = 0.101)") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB"))
  #facet_wrap(~Site, labeller = labeller(Site = labels_cover))

# Save plot as a pdf
ggsave("NMDS_bact_cover.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# ITS

cover_fil_f_asv <- subset.data.frame(fil_f_asv[c("IOB1","IOB2", "IOB3", "IOB4", "IOB5",
                                                 "ION1", "ION2", "ION3", "ION4", "ION5",
                                                 "IOD1", "IOD2", "IOD3", "IOD4", "IOD5",
                                                 "NON1", "NON2", "NON3", "NON4", "NON5",
                                                 "NOD1", "NOD2", "NOD3", "NOD4", "NOD5", 
                                                 "COB1", "COB2", "COB3", "COB4", "COB5"),])


cover_f_asv_bray <- vegdist(cover_fil_f_asv, method="bray")
cover_f_asv_nmds <- metaMDS(cover_f_asv_bray, k=2, try = 100)
cover_f_metadata$Axis01 = cover_f_asv_nmds$points[,1]
cover_f_metadata$Axis02 = cover_f_asv_nmds$points[,2]

cover_f_asv_nmds$stress

ggplot(cover_f_metadata,aes(Axis01, Axis02))+
  geom_point(aes(color=Site, shape=SoilCover), size=4)+
  labs(title="ITS NMDS ordination")+
  theme_classic(base_size = 14)

# perform permanova
adonis2(cover_f_asv_bray ~ SoilCover*Site, data = cover_f_metadata)

# labels_cover <- c("I" = "Induced", "N" = "Natural Invasion")

# plot fungal NMDS
ggplot(data = cover_f_metadata, aes(Axis01, Axis02)) +
  geom_point(aes(color=SoilCover), size=4) +
  geom_mark_hull(aes(group = Site, label = Site), concavity = 10) +
  labs(title="Fungal NMDS ordination (Stress = 0.160)")+
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB"))
  #facet_wrap(~Site, labeller = labeller(Site = labels_cover))

# Save plot as a pdf
ggsave("NMDS_fung_cover.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# ---------------------------------------------------------------------------------
## Tree Analysis ##
# ---------------------------------------------------------------------------------

# data frame of only samples without trees
open_sample = c("IOB1","IOB2", "IOB3", "IOB4", "IOB5",
                "ION1", "ION2", "ION3", "ION4", "ION5",
                "IOD1", "IOD2", "IOD3", "IOD4", "IOD5",
                "NON1", "NON2", "NON3", "NON4", "NON5",
                "NOD1", "NOD2", "NOD3", "NOD4", "NOD5", 
                "COB1", "COB2", "COB3", "COB4", "COB5")

# subset of the bacteria metadata with only trees
trees_b_metadata <- b_metadata[!(row.names(b_metadata) %in% open_sample),]

# duplicate the vegetation column and call it "TreeSpecies"
trees_b_metadata$TreeSpecies = trees_b_metadata$Vegetation

# duplicate the Site column and call it "InvasionLevel"
trees_b_metadata$InvasionLevel = trees_b_metadata$Site


# subset of the fungi metadata with only trees
trees_f_metadata <- f_metadata[!(row.names(f_metadata) %in% open_sample),]

# duplicate the vegetation column and call it "TreeSpecies"
trees_f_metadata$TreeSpecies = trees_f_metadata$Vegetation

# duplicate the Site column and call it "InvasionLevel"
trees_f_metadata$InvasionLevel = trees_f_metadata$Site



# linear models & plots

# 16S

# plot of bacterial and archaeal richness (observed means and standard errors)
ggplot(data = trees_b_metadata, aes(x = InvasionLevel, y = Richness.16S.rar, color = TreeSpecies)) +
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data = "mean_se") +
  geom_jitter(aes(x=InvasionLevel, y=Richness.16S.rar), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Invasion Level") +
  ylab("Bacterial and Archaeal Richness") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_grid(.~factor(InvasionLevel, levels = c("N", "C", "I")), scales = "free", space = "free") +
  scale_x_discrete(labels=c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")) +
  theme(strip.text.x = element_blank())

# save plot as a pdf
ggsave("Rich_bact_trees.pdf", device = "pdf", width = 9, height = 7, units = "in")

# linear model with bacterial richness as response variable, and InvasionLevel and TreeSpecies as interacting explanatory variables 
trees_b_richness.lm <- lm(Richness.16S.rar ~ InvasionLevel * TreeSpecies, data = trees_b_metadata)
summary(trees_b_richness.lm)
anova(trees_b_richness.lm)

# estimate means of the previous model to plot them
trees_b_richness.means <- estimate_means(trees_b_richness.lm)
trees_b_richness.means

# plot of bacterial and archaeal richness (estimated means and confidence intervals)
ggplot(data = trees_b_richness.means, aes(x = InvasionLevel, y = Mean, color = TreeSpecies)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge2(width=1)) +
  geom_jitter(aes(x=InvasionLevel, y=Mean), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Invasion Level") +
  ylab("Bacterial and Archaeal Richness") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_grid(.~factor(InvasionLevel, levels = c("N", "C", "I")), scales = "free", space = "free") +
  scale_x_discrete(labels=c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")) +
  theme(strip.text.x = element_blank())

# Save plot as a pdf
ggsave("Rich_bact_trees_lm.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# plot of bacterial diversity (observed means and standard errors)
ggplot(data = trees_b_metadata, aes(x = InvasionLevel, y = Shannon.16S, color = TreeSpecies)) +
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data = "mean_se") +
  geom_jitter(aes(x=InvasionLevel, y=Shannon.16S), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Invasion Level") +
  ylab("Bacterial and Archaeal Diversity") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_grid(.~factor(InvasionLevel, levels = c("N", "C", "I")), scales = "free", space = "free") +
  scale_x_discrete(labels=c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")) +
  theme(strip.text.x = element_blank())

# save plot as a pdf
ggsave("Div_bact_trees.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# linear model with bacterial diversity as response variable
trees_b_diversity.lm <- lm(Shannon.16S ~ InvasionLevel * TreeSpecies, data = trees_b_metadata)
summary(trees_b_diversity.lm)
anova(trees_b_diversity.lm)

# estimate means of the previous model to plot them
trees_b_diversity.means <- estimate_means(trees_b_diversity.lm)
trees_b_diversity.means

# plot of bacterial and archaeal diversity (estimated means and confidence intervals)
ggplot(data = trees_b_diversity.means, aes(x = InvasionLevel, y = Mean, color = TreeSpecies)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge2(width=1)) +
  geom_jitter(aes(x=InvasionLevel, y=Mean), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Invasion Level") +
  ylab("Bacterial and Archaeal Diversity") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_grid(.~factor(InvasionLevel, levels = c("N", "C", "I")), scales = "free", space = "free") +
  scale_x_discrete(labels=c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")) +
  theme(strip.text.x = element_blank())

# Save plot as a pdf
ggsave("Div_bact_trees_lm.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# ITS

# plot of fungal richness (observed means and standard errors)
ggplot(data = trees_f_metadata, aes(x = InvasionLevel, y = Richness.ITS.rar, color = TreeSpecies)) +
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data = "mean_se") +
  geom_jitter(aes(x=InvasionLevel, y=Richness.ITS.rar), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Invasion Level") +
  ylab("Fungal Richness") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_grid(.~factor(InvasionLevel, levels = c("N", "C", "I")), scales = "free", space = "free") +
  scale_x_discrete(labels=c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")) +
  theme(strip.text.x = element_blank())

# save plot as pdf
ggsave("Rich_fung_trees.pdf", device = "pdf", width = 9, height = 7, units = "in")

# linear model with fungal richness as response variable
trees_f_richness.lm <- lm(Richness.ITS.rar ~ InvasionLevel * TreeSpecies, data = trees_f_metadata)
anova(trees_f_richness.lm)

# estimate means of the previous model to plot them
trees_f_richness.means <- estimate_means(trees_f_richness.lm)
trees_f_richness.means

# plot of fungal richness (estimated means and confidence intervals)
ggplot(data = trees_f_richness.means, aes(x = InvasionLevel, y = Mean, color = TreeSpecies)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge2(width=1)) +
  geom_jitter(aes(x=InvasionLevel, y=Mean), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Invasion Level") +
  ylab("Fungal Richness") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_grid(.~factor(InvasionLevel, levels = c("N", "C", "I")), scales = "free", space = "free") +
  scale_x_discrete(labels=c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")) +
  theme(strip.text.x = element_blank())

# Save plot as a pdf
ggsave("Rich_fung_trees_lm.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# plot of fungal diversity (observed means and standard errors)
ggplot(data = trees_f_metadata, aes(x = InvasionLevel, y = Shannon.ITS, color = TreeSpecies)) +
  geom_pointrange(position = position_dodge2(width=1), stat = "summary", fun.data = "mean_se") +
  geom_jitter(aes(x=InvasionLevel, y=Shannon.ITS), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Invasion Level") +
  ylab("Fungal Diversity") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_grid(.~factor(InvasionLevel, levels = c("N", "C", "I")), scales = "free", space = "free") +
  scale_x_discrete(labels=c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")) +
  theme(strip.text.x = element_blank())

# save plot as a pdf
ggsave("Div_fung_trees.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# linear model with fungal diversity as response variable
trees_f_diversity.lm <- lm(Shannon.ITS ~ InvasionLevel * TreeSpecies, data = trees_f_metadata)
summary(trees_f_diversity.lm)
anova(trees_f_diversity.lm)

# estimate means of the previous model to plot them
trees_f_diversity.means <- estimate_means(trees_f_diversity.lm)
trees_f_diversity.means

# plot of fungal diversity (estimated means and confidence intervals)
ggplot(data = trees_f_diversity.means, aes(x = InvasionLevel, y = Mean, color = TreeSpecies)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), position = position_dodge2(width=1)) +
  geom_jitter(aes(x=InvasionLevel, y=Mean), position=position_jitterdodge(jitter.width=0.1)) +
  xlab("Invasion Level") +
  ylab("Fungal Diversity") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_grid(.~factor(InvasionLevel, levels = c("N", "C", "I")), scales = "free", space = "free") +
  scale_x_discrete(labels=c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")) +
  theme(strip.text.x = element_blank())

# Save plot as a pdf
ggsave("Div_fung_trees_lm.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# Community composition 

# 16S

# data frame of only samples without trees
open_sample2 = c("iob1","iob2", "iob3", "iob4", "iob5",
                 "ion1", "ion2", "ion3", "ion4", "ion5",
                 "iod1", "iod2", "iod3", "iod4", "iod5",
                 "non1", "non2", "non3", "non4", "non5",
                 "nod1", "nod2", "nod3", "nod4", "nod5", 
                 "cob1", "cob2", "cob3", "cob4", "cob5")


trees_fil_b_asv <- fil_b_asv[!(row.names(fil_b_asv) %in% open_sample2),]


trees_b_asv_bray <- vegdist(trees_fil_b_asv, method="bray")
trees_b_asv_nmds <- metaMDS(trees_b_asv_bray, k=2, try = 100)
trees_b_metadata$Axis01 = trees_b_asv_nmds$points[,1]
trees_b_metadata$Axis02 = trees_b_asv_nmds$points[,2]

trees_b_asv_nmds$stress

ggplot(trees_b_metadata,aes(Axis01, Axis02))+
  geom_point(aes(color=InvasionLevel, shape=TreeSpecies), size=4)+
  labs(title="16S NMDS ordination")+
  theme_classic(base_size = 14)

# perform permanova
adonis2(trees_b_asv_bray ~ InvasionLevel*TreeSpecies, data = trees_b_metadata)

labels <- c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")

trees_b_metadata$InvasionLevel<-factor(trees_b_metadata$InvasionLevel,
                                       levels=c("N", "C", "I"))

# plot bacterial NMDS
ggplot(data = trees_b_metadata, aes(Axis01,Axis02)) +
  geom_point(aes(color=TreeSpecies), size=4) +
  labs(title="Bacterial NMDS ordination (Stress = 0.144)")+
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_wrap(~InvasionLevel, labeller = labeller(InvasionLevel = labels))

# save plot as a pdf
ggsave("NMDS_bact_trees.pdf", device = "pdf", width = 9, height = 7,
       units = "in")

# ITS

trees_fil_f_asv <- fil_f_asv[!(row.names(fil_f_asv) %in% open_sample2),]


trees_f_asv_bray <- vegdist(trees_fil_f_asv, method="bray")
trees_f_asv_nmds <- metaMDS(trees_f_asv_bray, k=2, try = 100)
trees_f_metadata$Axis01 = trees_f_asv_nmds$points[,1]
trees_f_metadata$Axis02 = trees_f_asv_nmds$points[,2]

trees_f_asv_nmds$stress

ggplot(trees_f_metadata,aes(Axis01, Axis02))+
  geom_point(aes(color=InvasionLevel, shape=TreeSpecies), size=4)+
  labs(title="ITS NMDS ordination")+
  theme_classic(base_size = 14)

# perform permanova
adonis2(trees_f_asv_bray ~ InvasionLevel*TreeSpecies, data = trees_f_metadata)

labels <- c("N" = "Non-invaded", "C" = "Invaded", "I" = "Induced")

trees_f_metadata$InvasionLevel<-factor(trees_f_metadata$InvasionLevel,
                                       levels=c("N", "C", "I"))

# plot fungal NMDS
ggplot(data = trees_f_metadata, aes(Axis01,Axis02)) +
  geom_point(aes(color=TreeSpecies), size=4) +
  labs(title="Fungal NMDS ordination (Stress = 0.236)") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41")) +
  facet_wrap(~InvasionLevel, labeller = labeller(InvasionLevel = labels))


# save plot as a pdf
ggsave("NMDS_fung_trees.pdf", device = "pdf", width = 9, height = 7 , units = "in")



# ---------------------------------------------------------------------------------
## Soil Data Analysis ##
# ---------------------------------------------------------------------------------

# PCA of soil data

# remove one of CPV2 (row 17)
soil_data.fil <- soil_data[-c(17),]

# change row names for the values under "Sample" column
row.names(soil_data.fil) <- soil_data.fil$Sample

# remove "Sample" column
soil_data.sub = subset(soil_data.fil, select = -c(Sample))

# PCA 
soil_data.pca <- rda(soil_data.sub, scale=T)
summary(soil_data.pca)

# extract eigenvalues
ev <- soil_data.pca$CA$eig

# screeplot
barplot(ev/sum(ev), main="Eigenvalues", col="bisque") 


plot(soil_data.pca, display="sites", xlab="PC1 (40%)", ylab="PC2 (27%)")
biplot(soil_data.pca, scaling="species", xlab="PC1 (40%)", ylab="PC2 (27%)") #Focus on variables (species)
biplot(soil_data.pca, scaling="sites", xlab="PC1 (40%)", ylab="PC2 (27%)") #Focus on objects (sites)

# Soil Cover #

# susbset of the soil metadata without trees
cover_soil_metadata <- subset.data.frame(soil_metadata.fil[c("IOB1","IOB2", "ION1", "ION2",
                                                             "IOD1", "IOD2","NON1", "NON2", "NOD1", 
                                                             "NOD2", "COB1", "COB2"),])

# replace site "C" with "N" in order to have 2 sites instead of 3
cover_soil_metadata$Site[cover_soil_metadata$Site == "C"] <- "N" 

# duplicate the vegetation column and call it "SoilCover"
cover_soil_metadata$SoilCover = cover_soil_metadata$Vegetation

# linear model with water content as response variable
cover_soil_water.lm <- lm(WaterContent ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_water.lm)
anova(cover_soil_water.lm)

# linear model with pH as response variable
cover_soil_pH.lm <- lm(pH ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_pH.lm)
anova(cover_soil_pH.lm)

# linear model with EC as response variable
cover_soil_EC.lm <- lm(EC ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_EC.lm)
anova(cover_soil_EC.lm)

# linear model with N % as response variable
cover_soil_N.lm <- lm(N ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_N.lm)
anova(cover_soil_N.lm)

# linear model with C % as response variable
cover_soil_C.lm <- lm(C ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_C.lm)
anova(cover_soil_C.lm)

# linear model with Fe as response variable
cover_soil_Fe.lm <- lm(Fe ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_Fe.lm)
anova(cover_soil_Fe.lm)

# linear model with Cu as response variable
cover_soil_Cu.lm <- lm(Cu ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_Cu.lm)
anova(cover_soil_Cu.lm)

# linear model with Zn as response variable
cover_soil_Zn.lm <- lm(Zn ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_Zn.lm)
anova(cover_soil_Zn.lm)

# linear model with K as response variable
cover_soil_K.lm <- lm(K ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_K.lm)
anova(cover_soil_K.lm)

# linear model with Mg as response variable
cover_soil_Mg.lm <- lm(Mg ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_Mg.lm)
anova(cover_soil_Mg.lm)

# linear model with Ca as response variable
cover_soil_Ca.lm <- lm(Ca ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_Ca.lm)
anova(cover_soil_Ca.lm)

# linear model with P as response variable
cover_soil_P.lm <- lm(P ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_P.lm)
anova(cover_soil_P.lm)

# linear model with Mn as response variable
cover_soil_Mn.lm <- lm(Mn ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_Mn.lm)
anova(cover_soil_Mn.lm)

# linear model with S as response variable
cover_soil_S.lm <- lm(S ~ Site * SoilCover, data = cover_soil_metadata)
summary(cover_soil_S.lm)
anova(cover_soil_S.lm)

# PCA of only soil cover 

# subset of soil data without trees
soil_data_cover <- subset.data.frame(soil_data.sub[c("IOB1","IOB2", "ION1", "ION2",
                                                             "IOD1", "IOD2","NON1", "NON2", "NOD1", 
                                                             "NOD2", "COB1", "COB2"),])

soil_data_cover.pca <- rda(soil_data_cover, scale=T)
summary(soil_data_cover.pca)

cover_soil_pca_metadata<-cbind(cover_soil_metadata, scores(soil_data_cover.pca)$sites)

ggplot(data = cover_soil_pca_metadata, aes(PC1, PC2)) +
  geom_point(aes(color=SoilCover), size=4) +
  geom_mark_hull(aes(group = Site, label = Site), concavity = 10) +
  labs(title="Soil Cover PCA") +
  scale_color_manual(values=c("#BFAB37", "#555599", "#66BBBB"))

# Save plot as a pdf
ggsave("PCA_soil_cover.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# Tree Analysis

# data frame of only samples without trees
open_sample_soil = c("IOB1","IOB2", "ION1", "ION2",
                     "IOD1", "IOD2","NON1", "NON2", "NOD1", 
                     "NOD2", "COB1", "COB2")

# subset of the soil metadata with only trees
trees_soil_metadata <- soil_metadata.fil[!(row.names(soil_metadata.fil) %in% open_sample_soil),]

# duplicate the vegetation column and call it "TreeSpecies"
trees_soil_metadata$TreeSpecies = trees_soil_metadata$Vegetation

# duplicate the Site column and call it "InvasionLevel"
trees_soil_metadata$InvasionLevel = trees_soil_metadata$Site


# linear model with water content as response variable
trees_soil_water.lm <- lm(WaterContent ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_water.lm)
anova(trees_soil_water.lm)

# linear model with pH as response variable
trees_soil_pH.lm <- lm(pH ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_pH.lm)
anova(trees_soil_pH.lm)

# linear model with EC as response variable
trees_soil_EC.lm <- lm(EC ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_EC.lm)
anova(trees_soil_EC.lm)

# linear model with N % as response variable
trees_soil_N.lm <- lm(N ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_N.lm)
anova(trees_soil_N.lm)

# linear model with C % as response variable
trees_soil_C.lm <- lm(C ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_C.lm)
anova(trees_soil_C.lm)

# linear model with Fe as response variable
trees_soil_Fe.lm <- lm(Fe ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_Fe.lm)
anova(trees_soil_Fe.lm)

# linear model with Cu as response variable
trees_soil_Cu.lm <- lm(Cu ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_Cu.lm)
anova(trees_soil_Cu.lm)

# linear model with Zn as response variable
trees_soil_Zn.lm <- lm(Zn ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_Zn.lm)
anova(trees_soil_Zn.lm)

# linear model with K as response variable
trees_soil_K.lm <- lm(K ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_K.lm)
anova(trees_soil_K.lm)

# linear model with Mg as response variable
trees_soil_Mg.lm <- lm(Mg ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_Mg.lm)
anova(trees_soil_Mg.lm)

# linear model with Ca as response variable
trees_soil_Ca.lm <- lm(Ca ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_Ca.lm)
anova(trees_soil_Ca.lm)

# linear model with P as response variable
trees_soil_P.lm <- lm(P ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_P.lm)
anova(trees_soil_P.lm)

# linear model with Mn as response variable
trees_soil_Mn.lm <- lm(Mn ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_Mn.lm)
anova(trees_soil_Mn.lm)

# linear model with S as response variable
trees_soil_S.lm <- lm(S ~ InvasionLevel * TreeSpecies, data = trees_soil_metadata)
summary(trees_soil_S.lm)
anova(trees_soil_S.lm)

# PCA of only trees

# subset of soil data with trees
soil_data_trees <- soil_data.sub[!(row.names(soil_data.sub) %in% open_sample_soil),]

# PCA
soil_data_trees.pca <- rda(soil_data_trees, scale=T)
summary(soil_data_trees.pca)


trees_soil_pca_metadata<-cbind(trees_soil_metadata, scores(soil_data_trees.pca)$sites)

ggplot(data = trees_soil_pca_metadata, aes(PC1,PC2)) +
  geom_point(aes(color=TreeSpecies, shape=InvasionLevel), size=4) +
  labs(title="Trees PCA") +
  scale_color_manual(values=c("#729DA8", "#876F3B", "#A0BF41"))

# Save plot as a pdf
ggsave("PCA_soil_trees.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# ---------------------------------------------------------------------------------
## Soil Cover & Tree Analysis ##
# ---------------------------------------------------------------------------------


# BACTERIAL RICHNESS

# convert Site Invasion and CoverGeneral into factors
b_metadata$SiteInvasion<-as.factor(b_metadata$SiteInvasion)
b_metadata$CoverGeneral<-as.factor(b_metadata$CoverGeneral)
b_metadata$CoverGeneral <- factor(b_metadata$CoverGeneral, levels=c("Open", "Mesquite", "IronWood", "PaloVerde"))
b_metadata$Cover<-as.factor(b_metadata$Cover)
b_metadata$Cover <- factor(b_metadata$Cover, levels=c("Open_Bare", "Open_withNatives", "Open_withBuffel",
                                                      "Mesquite_noBuffel", "Mesquite_withBuffel",
                                                      "IronWood_noBuffel", "IronWood_withBuffel",
                                                      "PaloVerde_noBuffel", "PaloVerde_withBuffel"))

colors_cover<-c("Open_Bare"="burlywood", "Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

library(ggh4x) #for nested facet
ggplot(b_metadata, aes(x=Cover, y=Richness.16S.rar, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Bacterial richness") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

b_richness_lm<-lm(Richness.16S.rar ~ SiteInvasion * Cover, data=b_metadata)
anova(b_richness_lm)
library(modelbased)
estimate_contrasts(b_richness_lm, contrast = "SiteInvasion:Cover", p_adjust = "none")

# Save plot as a pdf
ggsave("Bacterial_Richness.pdf", device = "pdf", width = 12, height = 10, units = "in")

b_richness_lm<-lm(Richness.16S.rar ~ SiteInvasion * Cover, data = b_metadata)
anova(b_richness_lm)
library(modelbased)
estimate_contrasts(b_richness_lm, contrast = "SiteInvasion:Cover", p_adjust = "none")

# BACTERIAL COMPOSITION SIMILARITY
b_asv_bray <- vegdist(fil_b_asv, method="bray")
b_asv_nmds <- metaMDS(b_asv_bray, k=2, try = 100)
b_metadata$Axis01 = b_asv_nmds$points[,1]
b_metadata$Axis02 = b_asv_nmds$points[,2]

library(ggforce)
ggplot(b_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_hull(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10)+
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Bacterial_NMDS.pdf", device = "pdf", width = 12, height = 10, units = "in")

adonis2(b_asv_bray~SiteInvasion*Cover, data = b_metadata, permutations = 9999)

# Install pairwaiseAdonis
install.packages("devtools")
library(devtools)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis) 

b_metadata$SiteInvasion_Cover<-b_metadata$SiteInvasion:b_metadata$Cover #new variable for pairwise comparisons
pairwise.adonis2(b_asv_bray~SiteInvasion_Cover, data = b_metadata)


# FUNGAL RICHNESS


f_metadata$SiteInvasion<-as.factor(f_metadata$SiteInvasion)
f_metadata$CoverGeneral<-as.factor(f_metadata$CoverGeneral)
f_metadata$CoverGeneral <- factor(f_metadata$CoverGeneral, levels=c("Open", "Mesquite", "IronWood", "PaloVerde"))
f_metadata$Cover<-as.factor(f_metadata$Cover)
f_metadata$Cover <- factor(f_metadata$Cover, levels=c("Open_Bare", "Open_withNatives", "Open_withBuffel",
                                                      "Mesquite_noBuffel", "Mesquite_withBuffel",
                                                      "IronWood_noBuffel", "IronWood_withBuffel",
                                                      "PaloVerde_noBuffel", "PaloVerde_withBuffel"))

colors_cover<-c("Open_Bare"="burlywood", "Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

library(ggh4x) #for nested facet
ggplot(f_metadata, aes(x=Cover, y=Richness.ITS.rar, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Bacterial richness") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))


# Save plot as a pdf
ggsave("Fungal_Richness.pdf", device = "pdf",  width = 12, height = 10, units = "in")

f_richness_lm<-lm(Richness.ITS.rar ~ SiteInvasion * Cover, data = f_metadata)
anova(f_richness_lm)
library(modelbased)
estimate_contrasts(f_richness_lm, contrast = "SiteInvasion:Cover", p_adjust = "none")

# FUNGAL COMPOSITION SIMILARITY

f_asv_bray <- vegdist(fil_f_asv, method="bray")
f_asv_nmds <- metaMDS(f_asv_bray, k=2, try = 100)
f_metadata$Axis01 = f_asv_nmds$points[,1]
f_metadata$Axis02 = f_asv_nmds$points[,2]

library(ggforce)
ggplot(f_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_hull(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10)+
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Fungal_NMDS.pdf", device = "pdf", width = 12, height = 10, units = "in")

adonis2(f_asv_bray~SiteInvasion*Cover, data = f_metadata, permutations = 9999)

# Install pairwaiseAdonis
install.packages("devtools")
library(devtools)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis) 

f_metadata$SiteInvasion_Cover<-f_metadata$SiteInvasion:f_metadata$Cover #new variable for pairwise comparisons
pairwise.adonis2(f_asv_bray~SiteInvasion_Cover, data = f_metadata)

# SOIL PCA ORDINATION

# remove one of CPV2 (row 17)
soil_data.fil <- soil_data[-c(17),]

# change row names for the values under "Sample" column
row.names(soil_data.fil) <- soil_data.fil$Sample

# remove "Sample" column
soil_data.sub = subset(soil_data.fil, select = -c(Sample))

# PCA
library(ggcorrplot)
soil_corr<-cor(soil_data.sub, method="spearman")
soil_p.mat<-cor_pmat(soil_data.sub)

ggcorrplot(soil_corr, method = "circle", type = "upper", outline.color = "white",
                     colors = c("firebrick2", "white", "dodgerblue2"), p.mat = soil_p.mat, insig = "blank")


soil_data.pca <- rda(soil_data.sub, scale=T)
summary(soil_data.pca) #39.5% / 26.6%

soil_pca_metadata<-cbind(soil_metadata.fil, scores(soil_data.pca)$sites)
soil_pca_vectors<-data.frame(pc1 = scores(soil_data.pca)$species[,1], 
                             pc2 = scores(soil_data.pca)$species[,2])

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
