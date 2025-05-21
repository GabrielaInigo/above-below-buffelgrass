# Microbial data analysis for buffelgrass project 
# Author: Gabriela Iñigo
# Date: March 2023


# set working directory
setwd("/Users/gabrielainigo/Documents/GitHub/above-below-buffelgrass")


# LOAD DATA

# load metadata
metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)

# 16S 

# load bacterial taxonomic table
b_taxa = read.csv("buffelgrass_soilmicrobiome/data/16S/taxa_table_16S.csv", header = T)

# load bacterial ASV table
b_asv = read.csv("buffelgrass_soilmicrobiome/data/16S/asv_table_16S.csv", header = T)

# check out blanks
b_blanks = b_asv[c(1,2,3,4,5,6,7),]

# check if there is no contamination by adding the blanks
rowSums(b_blanks)

# remove blanks from asv table
b_asv = data.frame(b_asv[-c(1,2,3,4,5,6,7),])

# remove unwanted taxonomy from bacteria
install.packages("magrittr") # for %>% function
library(magrittr)
fil_b_taxa = b_taxa[!b_taxa$Order %in% "Chloroplast",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!fil_b_taxa$Family %in% "Mitochondria",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!fil_b_taxa$Kingdom %in% "Eukaryota",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!(is.na(fil_b_taxa$Kingdom)),] %>% droplevels()

# remove from ASV table
fil_b_asv <- b_asv[,rownames(fil_b_taxa)]

# change row names of filtered asv tables to upper case
# this is necessary to match rownames of asv tables with rownames of metadata
rownames(fil_b_asv) <- toupper(rownames(fil_b_asv))


# ITS

# load fungal taxonomic table 
f_taxa = read.csv("buffelgrass_soilmicrobiome/data/ITS/taxa_table_ITS.csv", header = T)

# load fungal ASV table
f_asv = read.csv("buffelgrass_soilmicrobiome/data/ITS/asv_table_ITS.csv", header = T)

# check out blanks ITS
f_blanks = f_asv[c(1,2,3,4,5,6,7),] 

# check if there is no contamination by adding the blanks
rowSums(f_blanks)

# remove blanks from asv table
f_asv = data.frame(f_asv[-c(1,2,3,4,5,6,7),])

# remove unwanted taxonomy from fungi
library(magrittr) # for %>% function
fil_f_taxa = f_taxa[f_taxa$Kingdom %in% "k__Fungi",] %>% droplevels() 

# remove from ASV table
fil_f_asv <- f_asv[,rownames(fil_f_taxa)]

# change row names of filtered asv tables to upper case
# this is necessary to match rownames of asv tables with rownames of metadata
rownames(fil_f_asv) <- toupper(rownames(fil_f_asv))


# DATA EXPLORATION

# 16S

# make histogram of bacterial ASV table to know the distribution of the data
hist(rowSums(fil_b_asv))

# get summary of the data
summary(rowSums(fil_b_asv)) # min = 40286

# get the standard deviation 
sd(rowSums(fil_b_asv)) # sd = 48106.9 

# order samples in ascending order
sort(rowSums(fil_b_asv)) # CM2 is the one with ~4000

# remove CM2 since it is the sample with the smallest number of ASVs
fil_b_asv  <- fil_b_asv[!(row.names(fil_b_asv) == "CM2"),]

# make a metadata for bacteria 
b_metadata <- metadata

# matches rownames metadata with rownames asv table
b_metadata = b_metadata[rownames(fil_b_asv),]

# make sure all are TRUE
rownames(b_metadata)==rownames(fil_b_asv) 

summary(rowSums(fil_b_asv)) # min = 40286

# create copy abundance matrix
fil_b_asv_ab <- fil_b_asv

# convert abundance to presence absence
# fil_b_asv_ab[fil_b_asv_ab > 0] <- 1

# calculate the mean of the total number of phylotypes in the samples
# mean(rowSums(fil_b_asv_ab)) # mean = 1664.473
# range(rowSums(fil_b_asv_ab)) # range = 612 - 2574
# summary(rowSums(fil_b_asv_ab)) # min = 612
# sd(rowSums(fil_b_asv_ab)) # sd = 490.5789

# ITS

# make histogram of fungal ASV table to know the distribution of the data
hist(rowSums(fil_f_asv)) 

# get summary of the data
summary(rowSums(fil_f_asv))

# get the standard deviation
sd(rowSums(fil_f_asv))

# order samples in ascending order
sort(rowSums(fil_f_asv)) # IOB5 is the one with ~1300

# remove iob5
fil_f_asv  <- fil_f_asv[!(row.names(fil_f_asv) == "IOB5"),]
f_metadata <- metadata
f_metadata = f_metadata[rownames(fil_f_asv),] #matches rownames metadata with rownames asv table
rownames(f_metadata)==rownames(fil_f_asv) #make sure all are TRUE

summary(rowSums(fil_f_asv))# min 27640

# create copy abundance matrix
fil_f_asv_ab <- fil_f_asv

# convert abundance to presence absence
# fil_f_asv_ab[fil_f_asv_ab > 0] <- 1

# calculate the mean of the total number of phylotypes in samples
# mean(rowSums(fil_f_asv_ab)) # mean = 267.1892
# range(rowSums(fil_f_asv_ab)) # range = 38 - 438
# summary(rowSums(fil_f_asv_ab)) # min = 38
# sd(rowSums(fil_f_asv_ab)) # sd = 100.0539


# TAXONOMY

# 16S

# install phyloseq package
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq", force = TRUE)

library(phyloseq)

# first we need to create a phyloseq objects of our data
# bacterial abundance table
count_tab_phy_b <- otu_table(as.matrix(fil_b_asv_ab), taxa_are_rows=F)

# bacterial taxa table
tax_tab_phy_b <- tax_table(as.matrix(fil_b_taxa))

# bacterial metadata table
sample_info_tab_phy_b <- sample_data(b_metadata)

b_ASV_physeq <- phyloseq(count_tab_phy_b, tax_tab_phy_b, sample_info_tab_phy_b)

#To get the ranking of classes: 
dat.aglo.b = tax_glom(b_ASV_physeq, taxrank = "Phylum")

# convert phyloseq object into a table
dat.dataframe.b = psmelt(dat.aglo.b) 

b.dat.agr = aggregate(Abundance~Phylum, data=dat.dataframe.b, FUN=sum)
b.dat.agr.sorted <- b.dat.agr[order(b.dat.agr$Abundance, decreasing=TRUE),]# this gives us the list of classes with their abundance

# calculate percentages of most abundant classes
b_proportions <- apply(b.dat.agr.sorted, 1, function(x){as.numeric(x)*100/sum(rowSums(b.dat.agr.sorted[2]))})

# Actinobacteria (35.97%)
# Proteobacteria (23.94%)
# Acidobacteria (9.22%) 
# Chloroflexi (8.28%)
# Firmicutes (6.83%)
# Gemmatimonadetes (3.20%)

sum(rowSums(b.dat.agr.sorted[2])) # 9778504

sum(b_proportions[2,]) # Must be 100%

# Taxonomic omposition barplot

# create phyloseq object
otu_b <- otu_table(as.matrix(fil_b_asv), taxa_are_rows = FALSE)
tax_b <- tax_table(as.matrix(fil_b_taxa))
sample_b <- sample_data(b_metadata)
b_physeq <- phyloseq(otu_b, tax_b, sample_b)

# agglomerate taxa
b_physeq_phylum <- tax_glom(b_physeq, taxrank = "Phylum")

# transform to relative abundance
b_physeq_phylum_rel <- transform_sample_counts(b_physeq_phylum, function(x) x / sum(x))

# plot taxonomic composition per sample
library(ggplot2)
plot_bar(b_physeq_phylum_rel, fill = "Phylum") +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# create a dataframe 
df_bac <- psmelt(b_physeq_phylum_rel)

# group less abundant phyla into a group called 'Others' (top 10)
top_n <- 6
phylum_sums <- df_bac %>% group_by(Phylum) %>% summarise(total = sum(Abundance)) %>% arrange(desc(total))
top_phyla <- phylum_sums$Phylum[1:top_n]
df_bac$Phylum_grouped <- ifelse(df_bac$Phylum %in% top_phyla, as.character(df_bac$Phylum), 'Others')

# Recalcular abundancia agrupada
df_bac_grouped <- df_bac %>% group_by(Sample, get(group_var), Phylum_grouped) %>% summarise(Abundance = sum(Abundance))
colnames(df_bac_grouped)[2] <- group_var


# ITS

# first we need to create a phyloseq objects of our data
# fungal abundance table
count_tab_phy_f <- otu_table(as.matrix(fil_f_asv_ab), taxa_are_rows=F)

# fungal taxa table
tax_tab_phy_f <- tax_table(as.matrix(fil_f_taxa))

# fungal metadat table
sample_info_tab_phy_f <- sample_data(f_metadata)

f_ASV_physeq <- phyloseq(count_tab_phy_f, tax_tab_phy_f, sample_info_tab_phy_f)

# to get the ranking of classes: 
dat.aglo.f = tax_glom(f_ASV_physeq, taxrank = "Phylum")

dat.dataframe.f = psmelt(dat.aglo.f) #melt it

f.dat.agr = aggregate(Abundance~Phylum, data=dat.dataframe.f, FUN=sum)
f.dat.agr.sorted <- f.dat.agr[order(f.dat.agr$Abundance, decreasing=TRUE),]# this gives us the list of classes with their abundance

# calculate percentages of most abundant classes
f_proportions <- apply(f.dat.agr.sorted, 1, function(x){as.numeric(x)*100/sum(rowSums(f.dat.agr.sorted[2]))})

# Ascomycota (79.57%) 
# Basidiomycota (12.40%) 
# Mortierellomycota (6.61%) 
# Glomeromycota (0.51%) 
# Chytridiomycota (0.44%)  
# Rozellomycota (0.37%)

sum(rowSums(f.dat.agr.sorted[2])) # 5092591

sum(f_proportions[2,]) # Must be 100

# ALPHA DIVERSITY

# 16S

# rarefy richness data to remove small outliers
library(vegan) # for specnumber function (vegan 2.6-4)
b_metadata$Richness.16S.rar <- specnumber(rrarefy(fil_b_asv, sample=40000)) 

sum(b_metadata$Richness.16S.rar) # 119,899
summary(b_metadata$Richness.16S.rar) # min = 608, max = 2,450

# rarefy diversity data to remove small outliers
b_metadata$Shannon.16S.rar <- diversity(rrarefy(fil_b_asv, sample=40000), index = "shannon")

# diversity data without rarefaction
b_metadata$Shannon.16S <- diversity(fil_b_asv)

# ITS

# rarify richness data to remove small outliers
library(vegan) # for specnumber function (vegan 2.6-4)
f_metadata$Richness.ITS.rar<-specnumber(rrarefy(fil_f_asv, sample=20000))

sum(f_metadata$Richness.ITS.rar) # 18,828
summary(f_metadata$Richness.ITS.rar) # min = 38, max = 404

# rarify diversity data to remove small outliers
f_metadata$Shannon.ITS.rar <- diversity(rrarefy(fil_f_asv, sample=20000), index = "shannon")

# diversity data without rarefaction
f_metadata$Shannon.ITS <- diversity(fil_f_asv)

# RICHNESS PLOTS

# 16S

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

# Bacterial richness plot
library(ggplot2) #for ggplot function
library(ggh4x) #for nested facet
ggplot(b_metadata, aes(x=Cover, y=Richness.16S.rar, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Bacterial/archaeal richness") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

# save plot as a pdf
ggsave("Bacterial_Richness.pdf", device = "pdf", width = 12, height = 10, units = "in")

# bacterial richness linear model
b_richness_lm<-lm(Richness.16S.rar ~ SiteInvasion * Cover, data = b_metadata)
anova(b_richness_lm)

library(performance)
r2(b_richness_lm)

# POST-HOC TEST 
#library(emmeans)
#emmeans(b_richness_lm,specs = Richness.16S.rar ~ SiteInvasion * Cover)

# Pairwise
library(modelbased) #for estimate_contrast function
library(dplyr)
estimates <- estimate_contrasts(b_richness_lm, contrast = "SiteInvasion:Cover", p_adjust = "none")
estimates <- as.data.frame(estimates)
estimates %>%
  filter(Level1 == "Induced Open_Bare",
         Level2 == "Induced Open_withNatives")
write.csv(estimates, file = "/Users/gabrielainigo/Documents/GitHub/bac-richness-est.csv",
          row.names = F)

# bacterial diversity plot
library(ggplot2) #for ggplot function
install.packages("ggh4x")
library(ggh4x) #for nested facet
ggplot(b_metadata, aes(x=Cover, y=Shannon.16S.rar, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Bacterial/archaeal diversity") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

# Save plot as a pdf
ggsave("Bacterial_Diversity.pdf", device = "pdf", width = 12, height = 10, units = "in")

# bacterial diversity linear model
b_diversity_lm<-lm(Shannon.16S.rar ~ SiteInvasion * Cover, data = b_metadata)
anova(b_diversity_lm)
library(performance)
r2(b_diversity_lm)

# Bacterial diversity plot without rarefaction
library(ggplot2) #for ggplot function
library(ggh4x) #for nested facet
ggplot(b_metadata, aes(x=Cover, y=Shannon.16S, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Bacterial/archaeal diversity") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

# Save plot as a pdf
ggsave("Bacterial_Diversity_NoRar.pdf", device = "pdf", width = 12, height = 10, units = "in")

# bacterial diversity without rarefaction linear model
b_diversity_norar_lm<-lm(Shannon.16S ~ SiteInvasion * Cover, data = b_metadata)
anova(b_diversity_norar_lm)
r2(b_diversity_norar_lm)

# ITS

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

# Fungal richness plot
library(ggplot2) #for ggplot function
library(ggh4x) #for nested facet
ggplot(f_metadata, aes(x=Cover, y=Richness.ITS.rar, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Fungal richness") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))


# Save plot as a pdf
ggsave("Fungal_Richness.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# Fungal richness linear model
f_richness_lm<-lm(Richness.ITS.rar ~ SiteInvasion * Cover, data = f_metadata)
anova(f_richness_lm)
r2(f_richness_lm)

library(modelbased) # for estimate_contrast fucntion
estimate_contrasts(f_richness_lm, contrast = "SiteInvasion:Cover", p_adjust = "none")

# Fungal diversity plot
library(ggplot2) #for ggplot function
library(ggh4x) #for nested facet
ggplot(f_metadata, aes(x=Cover, y=Shannon.ITS.rar, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Fungal diversity") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

# Save plot as a pdf
ggsave("Fungal_diversity.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# Fungal diversity linear model
f_diversity_lm<-lm(Shannon.ITS.rar ~ SiteInvasion * Cover, data = f_metadata)
anova(f_diversity_lm)
r2(f_diversity_lm)

# Fungal diversity plot without rarefaction
library(ggplot2) #for ggplot function
library(ggh4x) #for nested facet
ggplot(f_metadata, aes(x=Cover, y=Shannon.ITS, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Fungal diversity") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

# Save plot as a pdf
ggsave("Fungal_diversity_NoRar.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# Fungal diversity without rarefaction linear model
f_diversity_norar_lm<-lm(Shannon.ITS ~ SiteInvasion * Cover, data = f_metadata)
anova(f_diversity_norar_lm)
library(performance)
r2(f_diversity_norar_lm)

# COMPOSITION SIMILARITY NMDS PLOTS

# 16S

b_asv_bray <- vegdist(fil_b_asv, method="bray")
b_asv_nmds <- metaMDS(b_asv_bray, k=2, try = 100)
b_metadata$Axis01 = b_asv_nmds$points[,1]
b_metadata$Axis02 = b_asv_nmds$points[,2]

b_asv_nmds$stress # 0.1277071

library(ggforce)
ggplot(b_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_hull(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10)+
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Bacterial_NMDS.pdf", device = "pdf", width = 12, height = 10, units = "in")


# PERMANOVA
adonis2(b_asv_bray~SiteInvasion*Cover, data = b_metadata, permutations = 9999)

# Install pairwaiseAdonis
install.packages("devtools")
library(devtools)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

b_metadata$SiteInvasion_Cover<-b_metadata$SiteInvasion:b_metadata$Cover #new variable for pairwise comparisons
pairwise.adonis2(b_asv_bray~SiteInvasion_Cover, data = b_metadata)


# 16S NMDS with ellipses
library(vegan)

b_asv_bray <- vegdist(fil_b_asv, method="bray")
b_asv_nmds <- metaMDS(b_asv_bray, k=2, try = 100)
b_metadata$Axis01 = b_asv_nmds$points[,1]
b_metadata$Axis02 = b_asv_nmds$points[,2]

b_asv_nmds$stress # 0.1277071

library(ggforce)
ggplot(b_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_polygon(stat = "ellipse", aes(group=SiteInvasion, color="black"), fill=NA) +
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Bacterial_NMDS_ellipses.pdf", device = "pdf", width = 12, height = 10, units = "in")

# 16S NMDS with ellipses and labels
library(ggforce)
ggplot(b_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_ellipse(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10) +
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Bacterial_NMDS_ellipses.pdf", device = "pdf", width = 12, height = 10, units = "in")

# 16S NMDS with 95% confidence interval ellipses 
library(ggforce)
ggplot(b_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  stat_ellipse(geom = "polygon", aes(group=SiteInvasion, color="black"), fill=NA) +
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

ggsave("Bacterial_NMDS_ellipses_95.pdf", device = "pdf", width = 12, height = 10, units = "in")

# ITS

f_asv_bray <- vegdist(fil_f_asv, method="bray")
f_asv_nmds <- metaMDS(f_asv_bray, k=2, try = 100)
f_metadata$Axis01 = f_asv_nmds$points[,1]
f_metadata$Axis02 = f_asv_nmds$points[,2]

f_asv_nmds$stress # 0.1791455

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

#ITS NMDS with ellipses
f_asv_bray <- vegdist(fil_f_asv, method="bray")
f_asv_nmds <- metaMDS(f_asv_bray, k=2, try = 100)
f_metadata$Axis01 = f_asv_nmds$points[,1]
f_metadata$Axis02 = f_asv_nmds$points[,2]

f_asv_nmds$stress # 0.1791455

library(ggforce)
ggplot(f_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_polygon(stat = "ellipse", aes(group=SiteInvasion, color="black"), fill=NA) +
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Fungal_NMDS_ellipses.pdf", device = "pdf", width = 12, height = 10, units = "in")


#ITS NMDS with ellipses and labels
library(ggforce)
ggplot(f_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  geom_mark_ellipse(aes(group=SiteInvasion, color="black", label = SiteInvasion), concavity=10) +
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)
  
ggsave("Fungal_NMDS_ellipses.pdf", device = "pdf", width = 12, height = 10, units = "in")

#ITS NMDS with 95% confidence interval ellipses 
library(ggforce)
ggplot(f_metadata, aes(x=Axis01, y=Axis02, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4) +
  stat_ellipse(geom = "polygon", aes(group=SiteInvasion, color="black"), fill=NA) +
  scale_color_manual(values=colors_cover)+
  theme_classic(base_size = 14)

ggsave("Fungal_NMDS_ellipses_95.pdf", device = "pdf", width = 12, height = 10, units = "in")
