# preparing data for faprotax and fungild
# Author: Mery Touceda-Suárez
# Date: August 11th 2022

library(stringr)
library(tidyr)
library(see)
library(ggplot2)


setwd("/Volumes/BunnyBike/mex_buffel")

bac_asv <- read.csv("16S/asv_table_16S.csv", header = T)
bac_tax <- read.csv("16S/taxa_table_16S.csv", header = T)

fun_asv <- read.csv("ITS/asv_table_ITS.csv", header = T)
fun_tax <- read.csv("ITS/taxa_table_ITS.csv", header = T)

#-------------------------------------------------------------- Faprotax
# 1. Transpose asv table
bac_asv <- as.data.frame(t(bac_asv))

# 2. Join the columns in the taxa table into vectors that of words separated by |

bac_tax[is.na(bac_tax)]<- 'NA' #change nas for white spaces

bac_tax$combined <- str_c(bac_tax$Kingdom , ';', bac_tax$Phylum , ';', bac_tax$Class, ';', 
                          bac_tax$Order, ';', bac_tax$Family ,';', bac_tax$Genus, ';', 
                          bac_tax$Species)

# 3. Add taxonomy to ASV table 
bac_asv$taxonomy = bac_tax$combined[match(rownames(bac_asv),rownames(bac_tax))]

# 4. Write table

#bac_asv$ASV <- rownames(bac_asv)

write.table(bac_asv, '16S/asvtax_table_16S.txt', quote = F, sep = '\t')
fapro_input_table <- read.table('16S/asvtax_table_16S.txt', header = T)

fapro_input_table[,1]


# ------------------------------------------------------------- Fapro plotting!

fapro_filtered <- read.csv('16S/faprotax_table_filtered.csv', header = T,sep = ';')
rownames(fapro_filtered) <- fapro_filtered$group
fapro_filtered = fapro_filtered[,-1]

fapro_filtered <- as.data.frame(t(fapro_filtered))
bac_asv_t <- as.data.frame(t(bac_asv))
bac_asv_t <- bac_asv_t[,-c(1:7)]
sums = colSums(bac_asv_t)
fapro_filtered  = sweep(fapro_filtered , 1, sums, `/`)

metadata <- read.csv("mexbuff_metadata.csv")

test <- toupper(rownames(fapro_filtered))
rownames(fapro_filtered) <- test
metadata <- metadata[order(metadata$SampleID),]
metadata <- cbind(metadata, fapro_filtered)

metadata[metadata == "N"] <- "WithoutBuffel"
metadata[metadata == "C"] <- "WithBuffel"
metadata[metadata== "I"] <- "Induced"

levels(as.factor(metadata$Vegetation))

neworder <- c("Mesquite","Ironwood","Palo Verde", "Open Bare", "Open Buffel", "Open Native")
library(plyr)  ## or dplyr (transform -> mutate)
metadata <- arrange(transform(metadata,
                           Vegetation=factor(Vegetation,levels=neworder)),Vegetation)

ggplot(subset(metadata, Vegetation=="Mesquite" | Vegetation=="Ironwood" | Vegetation=="Palo Verde"), 
       aes(x=Site, y=photosynthetic_cyanobacteria, color=Site))+
 # geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  geom_pointrange(position=position_nudge(x=0), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of photosynthetic cyanobacteria") +
  scale_x_discrete(limits = rev(levels(as.factor(metadata$Site))))+
  scale_color_manual(values=c("firebrick4", "firebrick1", "lightblue3"))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  theme(legend.position="none", axis.text.x = element_text(size = 6))
ggsave("photosynthetic_cyanobacteria.pdf", device = "pdf", width = 9, height = 7 , units = "in")

colnames(metadata)




#-------------------------------------------------------------- Fungild

# 1. Transpose asv table
fun_asv <- as.data.frame(t(fun_asv))

# 2. Join the columns in the taxa table into vectors that of words separated by |

fun_tax[is.na(fun_tax)]<- 'NA' #change nas for white spaces

fun_tax$combined <- str_c(fun_tax$Kingdom , ';', fun_tax$Phylum , ';', fun_tax$Class, ';', 
              fun_tax$Order, ';', fun_tax$Family ,';', fun_tax$Genus, ';', 
              fun_tax$Species)

# 3. Add taxonomy to ASV table 
fun_asv$taxonomy = fun_tax$combined[match(rownames(fun_asv),rownames(fun_tax))]

# 4. Write table

#fun_asv$ASV <- rownames(fun_asv)

write.table(fun_asv, 'ITS/asvtax_table_ITS.txt', quote = F, sep = '\t')

library(FUNGuildR)


# --------------------------------------------- LET'S USE FUNGILD! 
#load previous table
fun_asv_tax <- read.table('ITS/asvtax_table_ITS.txt', header = T)
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
#calculate sum of abundances of asvs in each sample to make the ratio
sums = colSums(fun_asv[,1:82])
# make ratio
functional_table_final = sweep(functional_table_final, 2, sums, `/`)
functional_table_final = functional_table_final[c(4,10,12,22,27,35,36,41,44,49,52,54),]
functional_table_final = data.frame(functional_table_final[,8:82])

write.table(functional_table_final, "fungild_output_filtered.txt", quote = F, sep = ',')





# ------------------------------------------------------------ Plot Fungild!
function_table <- read.table("fungild_output_filtered.txt", header = T, sep = ',')
metadata <- read.csv("mexbuff_metadata.csv")

#merge function table with metadata
function_table <- as.data.frame(t(functional_table_final))
test <- toupper(rownames(function_table))
rownames(function_table) <- test
metadata <- metadata[order(metadata$SampleID),]
metadata <- cbind(metadata, function_table)

colnames(metadata)
#fix column names with spaces
colnames(metadata) <- sub(" ", "_", colnames(metadata))

metadata[metadata == "N"] <- "WithoutBuffel"
metadata[metadata == "C"] <- "WithBuffel"
metadata[metadata== "I"] <- "Induced"

neworder <- c("Mesquite","Ironwood","Palo Verde", "Open Bare", "Open Buffel", "Open Native")
library(plyr)  ## or dplyr (transform -> mutate)
metadata <- arrange(transform(metadata,
                              Vegetation=factor(Vegetation,levels=neworder)),Vegetation)

ggplot(subset(metadata, Vegetation=="Mesquite" | Vegetation=="Ironwood" 
              | Vegetation=="Palo Verde"), aes(x=Site, y=Wood_Saprotroph, color=Site))+
  #geom_pointrange(position=position_jitterdodge(), stat = "summary", fun.data="mean_se") +
  geom_pointrange(position=position_nudge(x=0), stat = "summary", fun.data="mean_se") +
  xlab(NULL) +
  ylab("Relative abundance of wood saprotrophs") +
  scale_x_discrete(limits = rev(levels(as.factor(metadata$Site))))+
  scale_color_manual(values=c("firebrick4", "firebrick1", "lightblue3"))+
  facet_wrap(~Vegetation)+
  theme_classic(base_size = 14)+
  scale_x_discrete(name = "",limits=c("WithoutBuffel", "WithBuffel", "Induced"))+
  theme(legend.position="none", axis.text.x = element_text(size = 6))
ggsave("Wood_Saprotroph.pdf", device = "pdf", width = 9, height = 7 , units = "in")

colnames(metadata)



