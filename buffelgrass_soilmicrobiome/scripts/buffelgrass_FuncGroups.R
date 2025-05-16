# Functional groups
# Author: Gabriela Iñigo 
# Date: May 2023


# set working directory
setwd("/Users/gabrielainigo/Documents/GitHub/")


# 16S 

# load bacterial ASV table
bac_asv <- read.csv("buffelgrass_soilmicrobiome/data/16S/asv_table_16S.csv", header = T)

# remove control blanks from bacterial ASV table
bac_asv = bac_asv[-c(1:7),]

# remove sample that was removed for bacterial richness
bac_asv <- bac_asv[!(row.names(bac_asv) == "cm2"),]

# change row names of bac_asv to uppercase
row.names(bac_asv) <- toupper(row.names(bac_asv))

# load filtered faprotax table
fapro_fil <- read.csv("buffelgrass_soilmicrobiome/data/16S/fapro_mexbuff_fil.csv", header = T)

# change rownmes of fapro_fil to the names of the functional groups in column X.groups
rownames(fapro_fil) <- fapro_fil$X.group 

# transpose fapro_fil
fapro_fil <- as.data.frame(t(fapro_fil))

# remove X.group column
fapro_fil = fapro_fil[-1,]

# change row names of fapro_fil to uppercase
row.names(fapro_fil) <- toupper(row.names(fapro_fil))

# make a new column called TotalSeqs in fapro_fil that has the total number of ASVs
fapro_fil$TotalSeqs <- rowSums(bac_asv)

# load metadata
metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)

# remove CM2 from metadata
metadata_b <- metadata[!(row.names(metadata) == "CM2"),]

# sort metadata_b row names in alphabetical order
metadata_b = metadata_b[order(row.names(metadata_b)), ]

# bind metadata with fapro_fil
fapro_metadata <- cbind(metadata_b, fapro_fil)

# make sure all are TRUE
rownames(fapro_metadata)==rownames(fapro_fil) 


# ITS

# load fungal ASV table
fun_asv <- read.csv("buffelgrass_soilmicrobiome/data/ITS/asv_table_ITS.csv", header = T)

# remove control blanks from fungal ASV table
fun_asv = fun_asv[-c(1:7),]

# remove sample that was removed for fungal richness
fun_asv <- fun_asv[!(row.names(fun_asv) == "iob5"),]

# change row names of fun_asv to uppercase
row.names(fun_asv) <- toupper(row.names(fun_asv))

# load filtered fungild table
fungild_fil <- read.csv("buffelgrass_soilmicrobiome/data/ITS/fungild_mexbuff_fil.csv", header = T)

# change rownmes of fungild_fil to the names of the functional groups in column Blank
rownames(fungild_fil) <- fungild_fil$Blank

# transpose fungild_fil
fungild_fil <- as.data.frame(t(fungild_fil))

# remove Blank row
fungild_fil = fungild_fil[-1,]

# remove X row
fungild_fil = fungild_fil[-75,]

# change row names of fungild_fil to uppercase
row.names(fungild_fil) <- toupper(row.names(fungild_fil))

# make a new column called TotalSeqs in fungild_fil that has the total number of ASVs
fungild_fil$TotalSeqs <- rowSums(fun_asv)

# make a new column called Saprotrophs in fungilf_fil that has the sum of all the different saprotrophs
fungild_fil$Saprotrophs <- as.numeric(fungild_fil$Dung_Saprotroph)+
  as.numeric(fungild_fil$Plant_Saprotroph)+
  as.numeric(fungild_fil$Undefined_Saprotroph)+
  as.numeric(fungild_fil$Wood_Saprotroph)

# load metadata
metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)

# remove IOB5 from metadata
metadata_f <- metadata[!(row.names(metadata) == "IOB5"),]

# sort metadata_f row names in alphabetical order
metadata_f = metadata_f[order(row.names(metadata_f)), ]

# bind metadata with fungild_fil
fungild_metadata <- cbind(metadata_f, fungild_fil)

# make sure all are TRUE
rownames(fungild_metadata)==rownames(fungild_fil) 



# Plot Bacterial Functional Groups

fapro_metadata$SiteInvasion<-as.factor(fapro_metadata$SiteInvasion)
fapro_metadata$CoverGeneral<-as.factor(fapro_metadata$CoverGeneral)
fapro_metadata$CoverGeneral <- factor(fapro_metadata$CoverGeneral, levels=c("Open", "Mesquite", "IronWood", "PaloVerde"))
fapro_metadata$Cover<-as.factor(fapro_metadata$Cover)
fapro_metadata$Cover <- factor(fapro_metadata$Cover, levels=c("Open_Bare", "Open_withNatives", "Open_withBuffel",
                                                              "Mesquite_noBuffel", "Mesquite_withBuffel",
                                                              "IronWood_noBuffel", "IronWood_withBuffel",
                                                              "PaloVerde_noBuffel", "PaloVerde_withBuffel"))

colors_cover<-c("Open_Bare"="burlywood", "Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

# Plot relative abundance of methylotrophs
library(ggplot2)
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y= as.numeric(methylotrophy)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Methylotrophs") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("methylotrophs.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change methylotrophy to numeric
fapro_metadata$methylotrophy <- as.numeric(fapro_metadata$methylotrophy)

# negative binomial model
library(MASS)
methyl_glm_nb <- glm.nb(methylotrophy ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(methyl_glm_nb)

library(performance) # for r2 function
r2(methyl_glm_nb)
anova(methyl_glm_nb, test = "Chisq")


# Plot relative abundance of nitrifiers
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(nitrification)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Nitrifiers") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("nitrifiers.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change nitrification to numeric
fapro_metadata$nitrification <- as.numeric(fapro_metadata$nitrification)


nitri_glm_p <- glm(nitrification ~ SiteInvasion * Cover, family = "poisson", data = fapro_metadata)
summary(nitri_glm_p)
r2(nitri_glm_p)
check_overdispersion(nitri_glm_p) #overdispersion

# since the data is overdispersed, a better model is a negative binomial glm 
library(MASS)
nitri_glm_nb <- glm.nb(nitrification ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(nitri_glm_nb)

r2(nitri_glm_nb)
anova(nitri_glm_nb, test = "Chisq")

# Plot relative abundance of denitrifiers
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(denitrification)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of denitrifiers") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("denitrifiers.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change denitrification to numeric
fapro_metadata$denitrification <- as.numeric(fapro_metadata$denitrification)

# negative binomial model
library(MASS)
denitri_glm_nb <- glm.nb(denitrification ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(denitri_glm_nb)

r2(denitri_glm_nb)
anova(denitri_glm_nb, test = "Chisq")


# Plot relative abundance of Chitinolytic Bacteria
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(chitinolysis)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Chitinolytic Bacteria") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("chitinolysis.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change chitinolysis to numeric
fapro_metadata$chitinolysis <- as.numeric(fapro_metadata$chitinolysis)

# negative binomial model
library(MASS)
chitin_glm_nb <- glm.nb(chitinolysis ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(chitin_glm_nb)

r2(chitin_glm_nb)
anova(chitin_glm_nb, test = "Chisq")

# Plot relative abundance of Nitrogen fixing bacteria
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(nitrogen_fixation)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative Abundance of Nitrogen Fixing Bacteria") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("nitrogen_fixation.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change nitrogen_fixation to numeric
fapro_metadata$nitrogen_fixation <- as.numeric(fapro_metadata$nitrogen_fixation)

# negative binomial model
library(MASS)
nitfix_glm_nb <- glm.nb(nitrogen_fixation ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(nitfix_glm_nb)

r2(nitfix_glm_nb)
anova(nitfix_glm_nb, test = "Chisq")

# Plot relative abundance of cellulolytic bacteria
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(cellulolysis)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of cellulolytic bacteria") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("cellulolysis.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change cellulolysis to numeric
fapro_metadata$cellulolysis <- as.numeric(fapro_metadata$cellulolysis)

# negative binomial model
library(MASS)
cellulo_glm_nb <- glm.nb(cellulolysis ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(cellulo_glm_nb)

r2(cellulo_glm_nb)
anova(cellulo_glm_nb, test = "Chisq")

# Plot relative abundance of xylanolytic bacteria
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(xylanolysis)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of xylanolytic bacteria") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("xylanolysis.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change xylanolysis to numeric
fapro_metadata$xylanolysis <- as.numeric(fapro_metadata$xylanolysis)

# negative binomial model
library(MASS)
xylano_glm_nb <- glm.nb(xylanolysis ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(xylano_glm_nb)

r2(xylano_glm_nb)
anova(xylano_glm_nb, test = "Chisq")

# Plot relative abundance of fermenting bacteria
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(fermentation)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of fermenting bacteria") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("fermentation.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change fermentation to numeric
fapro_metadata$fermentation <- as.numeric(fapro_metadata$fermentation)

# negative binomial model
library(MASS)
fermen_glm_nb <- glm.nb(fermentation ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(fermen_glm_nb)

r2(fermen_glm_nb)
anova(fermen_glm_nb, test = "Chisq")

# Plot relative abundance of aerobic chemoheterotrophic bacteria
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(aerobic_chemoheterotrophy)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of aerobic chemoheterotrophic bacteria") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("aerobic_chemoheterotrophy.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change aerobic_chemoheterotrophy to numeric
fapro_metadata$aerobic_chemoheterotrophy <- as.numeric(fapro_metadata$aerobic_chemoheterotrophy)

# negative binomial model
library(MASS)
aerochemohet_glm_nb <- glm.nb(aerobic_chemoheterotrophy ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(aerochemohet_glm_nb)

r2(aerochemohet_glm_nb)
anova(aerochemohet_glm_nb, test = "Chisq")

# Plot relative abundance of oxygenic photoautotrophic bacteria
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(oxygenic_photoautotrophy)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of oxygenic photoautotrophic bacteria") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("oxygenic_photoautotrophy.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change oxygenic_photoautotrophy to numeric
fapro_metadata$oxygenic_photoautotrophy <- as.numeric(fapro_metadata$oxygenic_photoautotrophy)

# negative binomial model
library(MASS)
oxyphotoaut_glm_nb <- glm.nb(oxygenic_photoautotrophy ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(oxyphotoaut_glm_nb)

r2(oxyphotoaut_glm_nb)
anova(oxyphotoaut_glm_nb, test = "Chisq")

# Plot relative abundance of ureolytic bacteria
library(ggh4x) #for nested facet
ggplot(fapro_metadata, aes(x=Cover, y=as.numeric(ureolysis)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of ureolytic bacteria") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("ureolysis.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change ureolysis to numeric
fapro_metadata$ureolysis <- as.numeric(fapro_metadata$ureolysis)

# negative binomial model
library(MASS)
ureol_glm_nb <- glm.nb(ureolysis ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fapro_metadata)
summary(ureol_glm_nb)

r2(ureol_glm_nb)
anova(ureol_glm_nb, test = "Chisq")


# Plot Fungal Functional Groups

fungild_metadata$SiteInvasion<-as.factor(fungild_metadata$SiteInvasion)
fungild_metadata$CoverGeneral<-as.factor(fungild_metadata$CoverGeneral)
fungild_metadata$CoverGeneral <- factor(fungild_metadata$CoverGeneral, levels=c("Open", "Mesquite", "IronWood", "PaloVerde"))
fungild_metadata$Cover<-as.factor(fungild_metadata$Cover)
fungild_metadata$Cover <- factor(fungild_metadata$Cover, levels=c("Open_Bare", "Open_withNatives", "Open_withBuffel",
                                                              "Mesquite_noBuffel", "Mesquite_withBuffel",
                                                              "IronWood_noBuffel", "IronWood_withBuffel",
                                                              "PaloVerde_noBuffel", "PaloVerde_withBuffel"))

colors_cover<-c("Open_Bare"="burlywood", "Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

# Plot relative abundance of Arbuscular Mycorrhizal Fungi
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Arbuscular_Mycorrhizal)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Arbuscular Mycorrhizal Fungi") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Arbuscular_Mycorrhizal.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Arbuscular_Mycorrhizal to numeric
fungild_metadata$Arbuscular_Mycorrhizal <- as.numeric(fungild_metadata$Arbuscular_Mycorrhizal)

# negative binomial model
library(MASS)
arbmyc_glm_nb <- glm.nb(Arbuscular_Mycorrhizal ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(arbmyc_glm_nb)

r2(arbmyc_glm_nb)
anova(arbmyc_glm_nb, test = "Chisq")

# Plot relative abundance of Dung Saprotroph Fungi
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Dung_Saprotroph)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Dung Saprotroph Fungi") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Dung_Saprotroph.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Dung_Saprotroph to numeric
fungild_metadata$Dung_Saprotroph <- as.numeric(fungild_metadata$Dung_Saprotroph)

# negative binomial model
library(MASS)
dungsapro_glm_nb <- glm.nb(Dung_Saprotroph ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(dungsapro_glm_nb)

r2(dungsapro_glm_nb)
anova(dungsapro_glm_nb, test = "Chisq")

# Plot relative abundance of Ectomycorrhizal Fungi
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Ectomycorrhizal)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Ectomycorrhizal Fungi") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Ectomycorrhizal.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Ectomycorrhizal to numeric
fungild_metadata$Ectomycorrhizal <- as.numeric(fungild_metadata$Ectomycorrhizal)

# negative binomial model
library(MASS)
ectomyco_glm_nb <- glm.nb(Ectomycorrhizal ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(ectomyco_glm_nb)

r2(ectomyco_glm_nb)
anova(ectomyco_glm_nb, test = "Chisq")

# Plot relative abundance of Lichenized Fungi
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Lichenized)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Lichenized Fungi") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Lichenized.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Lichenized to numeric
fungild_metadata$Lichenized <- as.numeric(fungild_metadata$Lichenized)

# negative binomial model
library(MASS)
lichen_glm_nb <- glm.nb(Lichenized ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(lichen_glm_nb)

r2(lichen_glm_nb)
anova(lichen_glm_nb, test = "Chisq")

# Plot relative abundance of Plant Pathogen Fungi
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Plant_Pathogen)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Plant Pathogen Fungi") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Plant_Pathogen.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Plant_Pathogen to numeric
fungild_metadata$Plant_Pathogen <- as.numeric(fungild_metadata$Plant_Pathogen)

# negative binomial model
library(MASS)
plantpath_glm_nb <- glm.nb(Plant_Pathogen ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(plantpath_glm_nb)

r2(plantpath_glm_nb)
anova(plantpath_glm_nb, test = "Chisq")

# Plot relative abundance of Plant Saprotroph Fungi
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Plant_Saprotroph)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Plant Saprotroph Fungi") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Plant_Saprotroph.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Plant_Saprotroph to numeric
fungild_metadata$Plant_Saprotroph <- as.numeric(fungild_metadata$Plant_Saprotroph)

# negative binomial model
library(MASS)
plantsapro_glm_nb <- glm.nb(Plant_Saprotroph ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(plantsapro_glm_nb)

r2(plantsapro_glm_nb)
anova(plantsapro_glm_nb, test = "Chisq")

# Plot relative abundance of Undefined Saprotroph Fungi
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Undefined_Saprotroph)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Undefined Saprotroph Fungi") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Undefined_Saprotroph.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Undefined_Saprotroph to numeric
fungild_metadata$Undefined_Saprotroph <- as.numeric(fungild_metadata$Undefined_Saprotroph)

# negative binomial model
library(MASS)
undefsapro_glm_nb <- glm.nb(Undefined_Saprotroph ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(undefsapro_glm_nb)

r2(undefsapro_glm_nb)
anova(undefsapro_glm_nb, test = "Chisq")

# Plot relative abundance of Wood Saprotroph Fungi
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Wood_Saprotroph)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Wood Saprotroph Fungi") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Wood_Saprotroph.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Wood_Saprotroph to numeric
fungild_metadata$Wood_Saprotroph <- as.numeric(fungild_metadata$Wood_Saprotroph)

# negative binomial model
library(MASS)
woodsapro_glm_nb <- glm.nb(Wood_Saprotroph ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(woodsapro_glm_nb)

r2(woodsapro_glm_nb)
anova(woodsapro_glm_nb, test = "Chisq")

# Plot relative abundance of Saprotrophs
library(ggh4x) #for nested facet
ggplot(fungild_metadata, aes(x=Cover, y= as.numeric(Saprotrophs)/TotalSeqs, color=Cover))+
  geom_jitter(alpha=0.6, position=position_jitterdodge(jitter.width = 0.2))+
  geom_pointrange(size=1, position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of Saprotrophs") +
  scale_color_manual(values=colors_cover)+
  facet_nested(~ SiteInvasion + CoverGeneral, scales = "free_x")+
  theme_classic(base_size = 14)+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# Save plot as a pdf
ggsave("Saprotrophs.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# change Saprotrophs to numeric
fungild_metadata$Saprotrophs <- as.numeric(fungild_metadata$Saprotrophs)

# negative binomial model
library(MASS)
sapro_glm_nb <- glm.nb(Saprotrophs ~ SiteInvasion * Cover + offset(log(TotalSeqs)), data = fungild_metadata)
summary(sapro_glm_nb)

r2(sapro_glm_nb)
anova(sapro_glm_nb, test = "Chisq")

#Scatterplots

#Scatterplot of nitrifiers and nitrogen fixers
library(ggforce)
ggplot(fapro_metadata, aes(x=as.numeric(nitrification)/TotalSeqs, y=as.numeric(nitrogen_fixation)/TotalSeqs, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4)+
  scale_color_manual(values=colors_cover)+
  xlab("Relative abundance of nitrifiers")+
  ylab("Relative abundance of nitrogen fixers")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Nitrifiers_Nfixers.pdf", device = "pdf", width = 12, height = 10, units = "in")

#Scatterplot of nitrifiers and denitrifiers
library(ggforce)
ggplot(fapro_metadata, aes(x=as.numeric(nitrification)/TotalSeqs, y=as.numeric(denitrification)/TotalSeqs, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4)+
  scale_color_manual(values=colors_cover)+
  xlab("Relative abundance of nitrifiers")+
  ylab("Relative abundance of denitrifiers")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Nitrifiers_Denitrifiers.pdf", device = "pdf", width = 12, height = 10, units = "in")

#Scatterplot of nitrifiers and ureolytic bacteria
library(ggforce)
ggplot(fapro_metadata, aes(x=as.numeric(nitrification)/TotalSeqs, y=as.numeric(ureolysis)/TotalSeqs, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4)+
  scale_color_manual(values=colors_cover)+
  xlab("Relative abundance of nitrifiers")+
  ylab("Relative abundance of ureolytic bacteria")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Nitrifiers_Ureolytic.pdf", device = "pdf", width = 12, height = 10, units = "in")

#Scatterplot of nitrogen fixers and denitrifiers
library(ggforce)
ggplot(fapro_metadata, aes(x=as.numeric(nitrogen_fixation)/TotalSeqs, y=as.numeric(denitrification)/TotalSeqs, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4)+
  scale_color_manual(values=colors_cover)+
  xlab("Relative abundance of nitrogen fixers")+
  ylab("Relative abundance of denitrifiers")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Nfixers_Denitrifiers.pdf", device = "pdf", width = 12, height = 10, units = "in")

#Scatterplot of nitrogen fixers and ureolytic bacteria
library(ggforce)
ggplot(fapro_metadata, aes(x=as.numeric(nitrogen_fixation)/TotalSeqs, y=as.numeric(ureolysis)/TotalSeqs, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4)+
  scale_color_manual(values=colors_cover)+
  xlab("Relative abundance of nitrogen fixers")+
  ylab("Relative abundance of ureolytic bacteria")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Nfixers_Ureolytic.pdf", device = "pdf", width = 12, height = 10, units = "in")

#Scatterplot of denitrifiers and ureolytic bacteria
library(ggforce)
ggplot(fapro_metadata, aes(x=as.numeric(denitrification)/TotalSeqs, y=as.numeric(ureolysis)/TotalSeqs, color=Cover))+
  geom_point(aes(shape=SiteInvasion), size=4)+
  scale_color_manual(values=colors_cover)+
  xlab("Relative abundance of denitrifiers")+
  ylab("Relative abundance of ureolytic bacteria")+
  theme_classic(base_size = 14)

# Save plot as a pdf
ggsave("Denitrifiers_Ureolytic.pdf", device = "pdf", width = 12, height = 10, units = "in")
