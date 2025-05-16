# Functional groups
# Author: Gabriela Iñigo 
# Date: May 2023


setwd("/Users/gabrielainigo/Documents/GitHub/")

bac_asv <- read.csv("buffelgrass_soilmicrobiome/data/16S/asv_table_16S.csv", header = T)
bac_tax <- read.csv("buffelgrass_soilmicrobiome/data/16S/taxa_table_16S.csv", header = T)

fun_asv <- read.csv("buffelgrass_soilmicrobiome/data/ITS/asv_table_ITS.csv", header = T)
fun_tax <- read.csv("buffelgrass_soilmicrobiome/data/ITS/taxa_table_ITS.csv", header = T)

# Bacterial Functional Groups (Faprotax)

# 1.Transpose asv table
bac_asv <- as.data.frame(t(bac_asv))


# 2.Join the columns in the taxa table into vectors that of words separated by |
bac_tax[is.na(bac_tax)]<- 'NA' #change NAs for white spaces

bac_tax$combined <- str_c(bac_tax$Kingdom , ';', bac_tax$Phylum , ';', bac_tax$Class, ';', 
                          bac_tax$Order, ';', bac_tax$Family ,';', bac_tax$Genus, ';', 
                          bac_tax$Species)


# 3.Add taxonomy to ASV table 
bac_asv$taxonomy = bac_tax$combined[match(rownames(bac_asv),rownames(bac_tax))]


# 4.Write table
write.table(bac_asv, 'buffelgrass_soilmicrobiome/data/16S/asvtax_table_16S.txt', quote = F, sep = '\t')
fapro_input_table <- read.table('buffelgrass_soilmicrobiome/data/16S/asvtax_table_16S.txt', header = T)

fapro_input_table[,1]

# Plot Bacterial Functional Groups

fapro_filtered <- read.csv('buffelgrass_soilmicrobiome/data/16S/fapro_mexbuff_fil.csv', header = T)

# change rownmes of fapro_filter to the names of the functional groups in column X.groups
rownames(fapro_filtered) <- fapro_filtered$X.group 

# remove X.group column
fapro_filtered = fapro_filtered[,-1]

# transpose fapro_filtered
fapro_filtered <- as.data.frame(t(fapro_filtered))

# make a new dataframe of bac_asv transpose
bac_asv_t <- as.data.frame(t(bac_asv))

# remove control blanks from bac_asv_t
bac_asv_t = bac_asv_t[-c(1:7),]


# transform all the values of the bac_asv_t data frame to numeric
bac_asv_t_num <- mutate_all(bac_asv_t, function(x) as.numeric(as.character(x)))

# make a variable with the sum of the rows of bac_asv_t (rows or columns?)
sums <- rowSums(bac_asv_t_num)

# get the relative abundance of each functional group by dividing it by the total sum of ASVs per sample (1 = rows)
fapro_filtered  = sweep(fapro_filtered, 1, sums, `/`)

# load metadata
metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)

# change row names to upper case for binding with metadata
test <- toupper(rownames(fapro_filtered))
rownames(fapro_filtered) <- test

# order the samples on the metadata for binding
metadata <- metadata[order(metadata$sampleid),]

# remove sample CM2 from metadata 
metadata = metadata[!(row.names(metadata) == "CM2"),]

# bind metadata and fapro_filtered
metadata_fapro <- cbind(metadata, fapro_filtered)

rownames(metadata_fapro)==rownames(fapro_filtered) #make sure all are TRUE



metadata_fapro$SiteInvasion<-as.factor(metadata_fapro$SiteInvasion)
metadata_fapro$CoverGeneral<-as.factor(metadata_fapro$CoverGeneral)
metadata_fapro$CoverGeneral <- factor(metadata_fapro$CoverGeneral, levels=c("Open", "Mesquite", "IronWood", "PaloVerde"))
metadata_fapro$Cover<-as.factor(metadata_fapro$Cover)
metadata_fapro$Cover <- factor(metadata_fapro$Cover, levels=c("Open_Bare", "Open_withNatives", "Open_withBuffel",
                                                      "Mesquite_noBuffel", "Mesquite_withBuffel",
                                                      "IronWood_noBuffel", "IronWood_withBuffel",
                                                      "PaloVerde_noBuffel", "PaloVerde_withBuffel"))

colors_cover<-c("Open_Bare"="burlywood", "Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

# Plot relative abundance of methylotrophs
library(ggh4x) #for nested facet
ggplot(metadata_fapro, aes(x=Cover, y=methylotrophy, color=Cover))+
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

# Plot relative abundance of nitrifiers
library(ggh4x) #for nested facet
ggplot(metadata_fapro, aes(x=Cover, y=nitrification, color=Cover))+
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


# Plot relative abundance of denitrifiers
library(ggh4x) #for nested facet
ggplot(metadata_fapro, aes(x=Cover, y=denitrification, color=Cover))+
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
ggsave("denitrifiers.pdf", device = "pdf",  width = 12, height = 10, units = "in")

# Plot relative abundance of Chitinolytic Bacteria
library(ggh4x) #for nested facet
ggplot(metadata_fapro, aes(x=Cover, y=chitinolysis, color=Cover))+
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


# Plot relative abundance of Nitrogen fixing bacteria
library(ggh4x) #for nested facet
ggplot(metadata_fapro, aes(x=Cover, y=nitrogen_fixation, color=Cover))+
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


# Plot relative abundance of cellulolytic bacteria
library(ggh4x) #for nested facet
ggplot(metadata_fapro, aes(x=Cover, y=cellulolysis, color=Cover))+
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

# Plot relative abundance of xylanolytic bacteria
library(ggh4x) #for nested facet
ggplot(metadata_fapro, aes(x=Cover, y=xylanolysis, color=Cover))+
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

# Plot relative abundance of fermenting bacteria
library(ggh4x) #for nested facet
ggplot(metadata_fapro, aes(x=Cover, y = fermentation, color=Cover))+
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



# Fungal Functional Groups (Fungild)


# 1.Transpose asv table
fun_asv <- as.data.frame(t(fun_asv))

# 2. Join the columns in the taxa table into vectors that of words separated by |
fun_tax[is.na(fun_tax)]<- 'NA' #change NAs for white spaces

fun_tax$combined <- str_c(fun_tax$Kingdom , ';', fun_tax$Phylum , ';', fun_tax$Class, ';', 
                          fun_tax$Order, ';', fun_tax$Family ,';', fun_tax$Genus, ';', 
                          fun_tax$Species)

# 3. Add taxonomy to ASV table 
fun_asv$taxonomy = fun_tax$combined[match(rownames(fun_asv),rownames(fun_tax))]

# 4. Write table
write.table(fun_asv, 'buffelgrass_soilmicrobiome/data/ITS/asvtax_table_ITS.txt', quote = F, sep = '\t')

# Install FUNGuildR 
devtools::install_github("brendanf/FUNGuildR")
library(FUNGuildR)

#load previous table
fun_asv_tax <- read.table('buffelgrass_soilmicrobiome/data/ITS/asvtax_table_ITS.txt', header = T)
#we only need the taxonomy column
#assign
assigned = funguild_assign(fun_asv_tax$taxonomy)
write.table(assigned, "fungild_assign_output.txt", quote = F)


#merge with ASV table
rownames(assigned) = rownames(fun_asv_tax)
fun_asv_functional = cbind(fun_asv, assigned)

###filter and keep only probable and highly probable
fun_asv_functional = fun_asv_functional[which(fun_asv_functional$confidenceRanking == "Probable"| 
                                                fun_asv_functional$confidenceRanking == "Highly Probable"), ]
####extract guilds
guilds = fun_asv_functional$guild


#take only asv abundances
functional_table_before_melt = data.frame(fun_asv_functional[,1:82])
#add the guilds and aggregate by sum of each guild
functional_table_before_melt$guilds = guilds
aggregated = aggregate(.~guilds, functional_table_before_melt, FUN = "sum")
#make final table with all but guild names
functional_table_final = aggregated[,-1]
#add guild names as rownames
rownames(functional_table_final)  = aggregated[,1]
write.table(functional_table_final, "fungild_output_aggregated.txt", quote = F, sep = ',')
#calculate sum of abundances of asvs in each sample to make the ratio
sums = colSums(fun_asv[,1:82])
# make ratio
functional_table_final = sweep(functional_table_final, 2, sums, `/`)
functional_table_final = functional_table_final[c(4,10,12,22,27,35,36,41,44,49,52,54),]
functional_table_final = data.frame(functional_table_final[,8:82])


write.table(functional_asv_final, "fungild_output_raw.txt", quote = F, sep = ',')
write.table(functional_table_final, "fungild_output_filtered.txt", quote = F, sep = ',')
