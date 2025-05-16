# Statistical analysis for the CAZ-MEX Buffelgrass soil microbiome project
# April 2022. 
# Author: Maria Touceda-Suarez

library(vegan)#multivar stats in ecology
library(data.table)
library(dbplyr) #data wrangling
library(tidyverse)
library(ggplot2) #pretty plots
library(lme4) #mixed effect models
library(lmerTest) #model testing for lme4
library(performance)
library(see)
library(ggforce)


setwd("/Users/gabrielainigo/Documents/GitHub/buffelgrass_soilmicrobiome/") #change directory to my path

#----------------------------------------------------------------------- 
#Load the data
#----------------------------------------------------------------------- 

#16S 
b_taxa = read.csv("data/16S/taxa_table_16S.csv", header=T)
#b_asv = fread("16S/asv_table_16S.csv", header=F, sep =',')
b_asv = read.csv("data/16S/asv_table_16S.csv", header=T)
metadata <- read.csv("data/mexbuff_metadata.csv", row.names = 1)

#ITS 
f_taxa = read.csv("data/ITS/taxa_table_ITS.csv", header=T)
#b_asv = fread("16S/asv_table_16S.csv", header=F, sep =',')
f_asv = read.csv("data/ITS/asv_table_ITS.csv", header=T)


# Check out blanks 16S
b_blanks = b_asv[c(1,2,3,4,5,6,7),] # numero de las columnas en b_asv que son controles
#check if they are 0
rowSums(b_blanks) # la suma de las rows que son los controles son para comprabar que no hay contaminacion
#  Blank   Zymo     b1     b2     b3     b4     b5 
# 2 167036     99     85     21     77   1076 

# Check out blanks ITs
f_blanks = f_asv[c(1,2,3,4,5,6,7),] 
#check if they are 0
rowSums(f_blanks)
# Blank  Zymo    b1    b2    b3    b4    b5 
# 14 47395 80601   959   738   764     0 

# remove blanks from asv table
b_asv = data.frame(b_asv[-c(1,2,3,4,5,6,7),])
f_asv = data.frame(f_asv[-c(1,2,3,4,5,6,7),])

# remove not wanted taxonomy
fil_b_taxa = b_taxa[!b_taxa$Order %in% "Chloroplast",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!fil_b_taxa$Family %in% "Mitochondria",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!fil_b_taxa$Kingdom %in% "Eukaryota",] %>% droplevels()
fil_b_taxa = fil_b_taxa[!(is.na(fil_b_taxa$Kingdom)),] %>% droplevels() #if kingdom is bact/arch, probably a real thing, although unknown taxonomy after that
# ! = - (funcionan igual)

#ITS
###sometimes there are "k__Eukaryota_kgd_Incertae_sedis" "k__Metazoa" "k__Rhizaria"
fil_f_taxa = f_taxa[f_taxa$Kingdom %in% "k__Fungi",] %>% droplevels() #puede amplificar otras cosa

# remove from ASV table
fil_b_asv<-b_asv[,rownames(fil_b_taxa)]
fil_f_asv<-f_asv[,rownames(fil_f_taxa)]


#----------------------------------------------------------------------- 
# Data exploration
#----------------------------------------------------------------------- 
# 1. Summary 
#16S
hist(rowSums(fil_b_asv)) # distribucion de la abundancia total de cada row
summary(rowSums(fil_b_asv))
sd(rowSums(fil_b_asv))

sort(rowSums(fil_b_asv)) # cm2 is the one with ~4000

# Remove cm2
fil_b_asv  <- fil_b_asv[!(row.names(fil_b_asv) == "cm2"),] # asv - amplicon sequence variant
b_metadata <- metadata
b_metadata = b_metadata[!metadata$sampleid %in% "cm2",] %>% droplevels()

summary(rowSums(fil_b_asv)) # min = 40286



#ITS
hist(rowSums(fil_f_asv)) 
summary(rowSums(fil_f_asv))
sd(rowSums(fil_f_asv))
sort(rowSums(fil_f_asv)) # iob5 is the one with ~1300

# Remove iob5
fil_f_asv  <- fil_f_asv[!(row.names(fil_f_asv) == "iob5"),]
f_metadata <- metadata
f_metadata = f_metadata[!metadata$sampleid %in% "iob5",] %>% droplevels()

summary(rowSums(fil_f_asv))# min 27640



# 2. Number of sequences (sum of reads over ASVs and samples)

# 3. Average number of phylotypes per sample
#16S
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



# ---------------------------------------------------------------------------------
## ALPHA DIVERSITY ##
# ---------------------------------------------------------------------------------

# 16S
b_metadata$Richness.16S.rar<-specnumber(rrarefy(fil_b_asv, sample=40000)) # remove small outliers
b_metadata$Shannon.16S <- diversity(rrarefy(fil_b_asv, sample=40000), index = "shannon")

# By site: I = (induced), C = (invasion + buffel), N = (invasion - buffel) 
ggplot(b_metadata, aes(x=Site, y=Richness.16S.rar)) +
  geom_jitter(position=position_jitterdodge(), aes(color=Site), size=2)+
  geom_boxplot(aes(fill=Site), outlier.shape = NA, alpha = 0.4)+
  xlab(NULL) +
  ylab("Bacterial/archaeal phylotype richness") +
  scale_x_discrete(labels= c(" (C) invasion + buffel"," (I) induced","(N) invasion - buffel"))+
  scale_fill_manual(values = c("tomato2", "orangered4", "goldenrod"))+ #for category
  scale_color_manual(values = c("tomato2", "orangered4", "goldenrod"))+ #for category
  theme_classic()+
  theme(legend.position="none",text = element_text(size=22))
#legend.position="none"
#ggsave("div_bact.pdf", device = "pdf", width = 9, height = 7 , units = "in")



ggplot(b_metadata, aes(x=Vegetation, y=Richness.16S.rar)) +
  geom_jitter(position=position_jitterdodge(), aes(color=Vegetation, shape = Site), size=2)+
  geom_boxplot(aes(fill=Vegetation), outlier.shape = NA, alpha = 0.4)+
  xlab(NULL) +
  ylab("Bacterial/archaeal phylotype richness") +
  #scale_x_discrete(limits=unique(md_sorted$Group))+
  scale_fill_manual(values = c("gray64", "tomato4", "goldenrod4", "forestgreen" ,"goldenrod", "darkolivegreen3"))+ #for category
  scale_color_manual(values = c("gray64", "tomato4", "goldenrod4", "forestgreen" ,"goldenrod", "darkolivegreen3"))+ #for category
  theme_classic()+
  theme(legend.position="right",text = element_text(size=22), axis.text.x = element_text(angle=45, hjust=1, size=14))



# ITS
f_metadata$Richness.16S.rar<-specnumber(rrarefy(fil_f_asv, sample=20000))
f_metadata$Shannon.16S <- diversity(rrarefy(fil_f_asv, sample=20000), index = "shannon")

# By site: I = (induced), C = (invasion + buffel), N = (invasion - buffel) 
ggplot(f_metadata, aes(x=Site, y=Richness.16S.rar)) +
  geom_jitter(position=position_jitterdodge(), aes(color=Site), size=2)+
  geom_boxplot(aes(fill=Site), outlier.shape = NA, alpha = 0.4)+
  xlab(NULL) +
  ylab("Fungal phylotype richness") +
  scale_x_discrete(labels= c(" (C) invasion + buffel"," (I) induced","(N) invasion - buffel"))+
  scale_fill_manual(values = c("tomato2", "orangered4", "goldenrod"))+ #for category
  scale_color_manual(values = c("tomato2", "orangered4", "goldenrod"))+ #for category
  theme_classic()+
  theme(legend.position="none",text = element_text(size=22))
#legend.position="none"
#ggsave("div_bact.pdf", device = "pdf", width = 9, height = 7 , units = "in")



ggplot(f_metadata, aes(x=Vegetation, y=Richness.16S.rar)) +
  geom_jitter(position=position_jitterdodge(), aes(color=Vegetation, shape = Site), size=2)+
  geom_boxplot(aes(fill=Vegetation), outlier.shape = NA, alpha = 0.4)+
  xlab(NULL) +
  ylab("Fungal phylotype richness") +
  #scale_x_discrete(limits=unique(md_sorted$Group))+
  scale_fill_manual(values = c("gray64", "tomato4", "goldenrod4", "forestgreen" ,"goldenrod", "darkolivegreen3"))+ #for category
  scale_color_manual(values = c("gray64", "tomato4", "goldenrod4", "forestgreen" ,"goldenrod", "darkolivegreen3"))+ #for category
  theme_classic()+
  theme(legend.position="right",text = element_text(size=22), axis.text.x = element_text(angle=45, hjust=1, size=14))



# ---------------------------------------------------------------------------------
## ALPHA DIVERSITY FOR POSTER ##
# ---------------------------------------------------------------------------------

# 16S
b_metadata$Richness.16S.rar<-specnumber(rrarefy(fil_b_asv, sample=40000))
b_metadata$Shannon.16S <- diversity(rrarefy(fil_b_asv, sample=40000), index = "shannon")

b_metadata$Site<-as.factor(b_metadata$Site)
levels(b_metadata$Site)<-c("WithBuffel", "Induced", "WithoutBuffel")
b_metadata$Site <- factor(b_metadata$Site, levels=c("WithoutBuffel", "WithBuffel", "Induced"))

b_metadata$Vegetation<-as.factor(b_metadata$Vegetation)
levels(b_metadata$Vegetation)<-c("Ironwood", "Mesquite", "OpenBare", "OpenBuffel", "OpenNative", "PaloVerde")
b_metadata$Vegetation <- factor(b_metadata$Vegetation, levels=c("OpenNative", "OpenBare", "OpenBuffel", "Mesquite", "Ironwood", "PaloVerde"))

ggplot(b_metadata, aes(x=Vegetation, y=Richness.16S.rar, fill=Site))+
  geom_jitter(position=position_jitterdodge())+
  geom_boxplot(outlier.shape = NA, alpha=0.4)+
  xlab(NULL) +
  ylab("16S richness") +
  theme_classic(base_size = 14)

ggplot(b_metadata, aes(x=Vegetation, y=Richness.16S.rar, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("16S richness") +
  theme_classic(base_size = 14)

ggplot(b_metadata, aes(x=Site, y=Richness.16S.rar, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("16S richness") +
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)

#plot for poster
pdf("bacteria_richness.pdf", width=10)
ggplot(subset(b_metadata, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="PaloVerde"), 
       aes(x=Site, y=Richness.16S.rar, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Bacterial richness") +
  scale_color_manual(values=c("lightblue3", "firebrick1", "firebrick4" ))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none")
dev.off()

anova(lm(Richness.16S.rar~Vegetation+Site, data=subset(b_metadata, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="PaloVerde")))

# ITS
f_metadata$Richness.ITS.rar<-specnumber(rrarefy(fil_f_asv, sample=25000))
f_metadata$Shannon.ITS <- diversity(rrarefy(fil_f_asv, sample=25000), index = "shannon")

f_metadata$Site<-as.factor(f_metadata$Site)
levels(f_metadata$Site)<-c("WithBuffel", "Induced", "WithoutBuffel")
f_metadata$Site <- factor(f_metadata$Site, levels=c("WithoutBuffel", "WithBuffel", "Induced"))

f_metadata$Vegetation<-as.factor(f_metadata$Vegetation)
levels(f_metadata$Vegetation)<-c("Ironwood", "Mesquite", "OpenBare", "OpenBuffel", "OpenNative", "PaloVerde")
f_metadata$Vegetation <- factor(f_metadata$Vegetation, levels=c("OpenNative", "OpenBare", "OpenBuffel", "Mesquite", "Ironwood", "PaloVerde"))

ggplot(f_metadata, aes(x=Vegetation, y=Richness.ITS.rar, fill=Site))+
  geom_jitter(position=position_jitterdodge())+
  geom_boxplot(outlier.shape = NA, alpha=0.4)+
  xlab(NULL) +
  ylab("ITS richness") +
  theme_classic(base_size = 14)

ggplot(f_metadata, aes(x=Vegetation, y=Richness.ITS.rar, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("ITS richness") +
  theme_classic(base_size = 14)

ggplot(f_metadata, aes(x=Site, y=Richness.ITS.rar, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("ITS richness") +
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)

#plot for poster
pdf("fungi_richness.pdf", width=10)
ggplot(subset(f_metadata, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="PaloVerde"), 
       aes(x=Site, y=Richness.ITS.rar, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Fungal richness") +
  scale_color_manual(values=c("lightblue3", "firebrick1", "firebrick4" ))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none")
dev.off()

anova(lm(Richness.ITS.rar~Vegetation*Site, data=subset(f_metadata, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="PaloVerde")))

# ---------------------------------------------------------------------------------
## COMMUNITY SIMILARITY ##
# ---------------------------------------------------------------------------------

# 16S
b_asv_bray <- vegdist(fil_b_asv, method="bray")
b_asv_nmds <- metaMDS(b_asv_bray, k=2, try = 100)
b_metadata$Axis01 = b_asv_nmds$points[,1]
b_metadata$Axis02 = b_asv_nmds$points[,2]

b_asv_nmds$stress

ggplot(b_metadata,aes(Axis01, Axis02))+
  geom_point(aes(color=Site, shape=Vegetation), size=4)+
  labs(title="16S NMDS ordination")+
  theme_classic(base_size = 14)

adonis(b_asv_bray~Vegetation*Site, data = b_metadata)

#plot for poster
pdf("bacteria_ordination.pdf", width=10)
ggplot(subset(b_metadata, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="PaloVerde"),
       aes(Axis01, Axis02))+
  geom_point(aes(color=Site), size=4)+
  labs(title="Bacterial NMDS ordination (Stress = 0.128)")+
  scale_color_manual(values=c("lightblue3", "firebrick1", "firebrick4" ))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none")
dev.off()

# ITS
f_asv_bray <- vegdist(fil_f_asv, method="bray")
f_asv_nmds <- metaMDS(f_asv_bray, k=2, try = 100)
f_metadata$Axis01 = f_asv_nmds$points[,1]
f_metadata$Axis02 = f_asv_nmds$points[,2]

f_asv_nmds$stress

ggplot(f_metadata, aes(Axis01, Axis02))+
  geom_point(aes(color=Site, shape=Vegetation), size=4)+
  labs(title="ITS NMDS ordination")+
  theme_classic(base_size = 14)

ggplot(subset(f_metadata, Vegetation=="Mesquite" | Vegetation=="PaloVerde"), aes(Axis01, Axis02))+
  geom_point(aes(color=Site, shape=Vegetation), size=4)+
  labs(title="ITS NMDS ordination")+
  theme_classic(base_size = 14)

adonis(f_asv_bray~Vegetation*Site, data = f_metadata)

#plot for poster
pdf("fungi_ordination.pdf", width=10)
ggplot(subset(f_metadata, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="PaloVerde"),
       aes(Axis01, Axis02))+
  geom_point(aes(color=Site), size=4)+
  labs(title="Fungal NMDS ordination (Stress = 0.179)")+
  scale_color_manual(values=c("lightblue3", "firebrick1", "firebrick4" ))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none")
dev.off()


# ---------------------------------------------------------------------------------
## VEGETATION ##
# ---------------------------------------------------------------------------------

veg <- read.csv("/Users/gabrielainigo/Documents/GitHub/buffelgrass_vegetation/VegAbundanceTable.csv", row.names = 1)

veg <- veg[,-1]

hist(rowSums(veg))
summary(rowSums(veg))
veg <- veg[order(rownames(metadata_veg)),]

metadata_veg <- metadata[!metadata$Vegetation %in% c("Open Buffel", "Open Native", "Open Bare"),]
# rarefaction?
metadata_veg$Rich.Veg <- specnumber(veg, MARGIN = 1)
metadata_veg$Shannon.Veg <- diversity(veg, MARGIN = 1)

library(ggplot2)
ggplot(subset(metadata_veg, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="Palo Verde"), 
       aes(x=Site, y=Rich.Veg, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Plant richness") +
  scale_color_manual(values=c("lightblue3", "firebrick1", "firebrick4" ))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none")

ggplot(subset(metadata_veg, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="Palo Verde"), 
       aes(x=Site, y=Shannon.Veg, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Plant Shannon diversity") +
  scale_color_manual(values=c("lightblue3", "firebrick1", "firebrick4" ))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none")

anova(lm(Rich.Veg~Vegetation+Site, data=subset(metadata_veg, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="Palo Verde")))
# que modelo usamos teniendo en cuenta los ceros? Podemos usar el mismo mdelo que para la microbial rich?



# asignar solo abundance buffel
# no encontrmos bufel!!!? NO LO INCLUYERON - si lo incluyeron pero como cencil (Cenchrus ciliare) y no pencil (Pennisetum ciliare)

metadata_veg$buffel_abundance <- veg$cencil
colSums(veg)
library(see)

ggplot(subset(metadata_veg, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="Palo Verde"), 
       aes(x=Site, y=buffel_abundance, color=Site))+
  geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Plant Shannon diversity") +
  scale_color_manual(values=c("lightblue3", "firebrick1", "firebrick4" ))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none")



# ordination
veg_bray <- vegdist(veg, method="bray")
veg_nmds <- metaMDS(veg_bray, k=2, try = 100)
metadata_veg$Axis01 = veg_nmds$points[,1]
metadata_veg$Axis02 = veg_nmds$points[,2]

veg_nmds$stress


ggplot(subset(metadata_veg, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="Palo Verde"),
       aes(Axis01, Axis02))+
  geom_point(aes(color=Site), size=4)+
  labs(title="Vegetation NMDS ordination (Stress = 0.179)")+
  scale_color_manual(values=c("lightblue3", "firebrick1", "firebrick4" ))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none")

adonis(veg_bray~Vegetation*Site, data = metadata_veg)


