# ITS NMDS Plot with correct package versions
# Author: Gabriela Iñigo
# Date: March 2024

# First, remove all packages from the search path
search()
lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""), detach, character.only=TRUE, unload=TRUE)

# Install and load required packages in correct order
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggforce", quietly = TRUE)) {
  install.packages("ggforce")
}
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}

# Load packages in correct order
library(ggplot2)
library(ggforce)
library(ggnewscale)

# Load the data
f_metadata <- read.csv("buffelgrass_soilmicrobiome/data/mexbuff_metadata.csv", row.names = 1)
fil_f_asv <- read.csv("buffelgrass_soilmicrobiome/data/ITS/asv_table_ITS.csv", header = T)

# Calculate NMDS
f_asv_bray <- vegdist(fil_f_asv, method="bray")
f_asv_nmds <- metaMDS(f_asv_bray, k=2, try = 100)
f_metadata$Axis01 = f_asv_nmds$points[,1]
f_metadata$Axis02 = f_asv_nmds$points[,2]

# Print stress value
f_asv_nmds$stress

# Set up factors for plotting
f_metadata$SiteInvasion <- as.factor(f_metadata$SiteInvasion)
f_metadata$Cover <- as.factor(f_metadata$Cover)

# Define colors for cover types
colors_cover <- c("Open_Bare"="burlywood", "Open_withNatives"="chocolate1", "Open_withBuffel"="chocolate4",
                 "Mesquite_noBuffel"="dodgerblue", "Mesquite_withBuffel"="dodgerblue4",
                 "IronWood_noBuffel"="firebrick1", "IronWood_withBuffel"="firebrick3",
                 "PaloVerde_noBuffel"="chartreuse2", "PaloVerde_withBuffel"="chartreuse4")

# Which Cover levels should be filled
filled_levels <- c("Mesquite_noBuffel", 
                   "IronWood_noBuffel", 
                   "PaloVerde_noBuffel", 
                   "Open_Bare", 
                   "Open_withNatives")

# Get all Cover levels present in your data
all_levels <- unique(f_metadata$Cover)

# Create a named vector for fill: matching the outline color if in filled_levels, otherwise NA
fill_values <- setNames(rep(NA, length(all_levels)), all_levels)
fill_values[filled_levels] <- colors_cover[filled_levels]

# Create the NMDS plot with filled points and ellipses
p <- ggplot(f_metadata, aes(x = Axis01, y = Axis02)) +
  # Points: color and fill both map to Cover
  geom_point(
    aes(color = Cover, fill = Cover, shape = SiteInvasion),
    size = 4
  ) +
  # First color/fill scales for Cover
  scale_color_manual(
    name   = "Cover",
    values = colors_cover
  ) +
  scale_fill_manual(
    name      = "Cover",
    values    = fill_values,
    na.value  = NA
  ) +
  # Reset color scale for ellipses
  new_scale_color() +
  # Ellipses: color by SiteInvasion
  geom_polygon(
    stat  = "ellipse",
    aes(group = SiteInvasion, color = SiteInvasion),
    fill  = NA,
    size  = 1
  ) +
  # Second color scale for SiteInvasion
  scale_color_manual(
    name   = "SiteInvasion",
    values = c("Induced" = "blue",
               "Natural" = "red")
  ) +
  # Control shapes by SiteInvasion
  scale_shape_discrete(name = "SiteInvasion") +
  # Overall theme
  theme_classic(base_size = 14) +
  # Add axis labels
  labs(x = "NMDS1", y = "NMDS2")

# Save plot as a pdf
ggsave("Fungal_NMDS_new.pdf", device = "pdf", width = 12, height = 10, units = "in") 