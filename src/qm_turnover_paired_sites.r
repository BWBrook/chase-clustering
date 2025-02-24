# Computes and plots site-specific compositional turnover for 34 sites that have both Pleistocene and Holocene data

import::from("dplyr", rowwise, mutate, select, filter)
import::from("magrittr", "%>%")
import::from("tibble", rownames_to_column)
library(ggplot2)
library(sf)

setwd("C:/git/chase-clustering") # adapt as needed

# Load the data
pleistocene <- read.csv("data/pleistocene_paired.csv", row.names = 1, check.names = FALSE)
holocene    <- read.csv("data/holocene_paired.csv",    row.names = 1, check.names = FALSE)
coords      <- read.csv("data/paired_site_coords.csv")

# Ensure both occupancy matrices have the same columns (sites)
# They should both have columns c_1, c_2, ..., c_34
pleistocene <- pleistocene[, coords$site, drop = FALSE]
holocene    <- holocene[, coords$site, drop = FALSE]

# Convert to presence/absence matrices in case there's any numerical weirdness
pleistocene <- (pleistocene > 0) * 1
holocene    <- (holocene > 0) * 1

# For each site (column), compute turnover between Pleistocene and Holocene
# A standard measure is "turnover" = (lost_species + gained_species) / (total_species across the two periods)
# turnover[i] = (# of species in Pleistocene only + # of species in Holocene only) / (# of species in Pleistocene union Holocene)
turnover_df <- data.frame(site = coords$site, turnover = NA_real_)

for (i in seq_along(coords$site)) {
  sitename <- coords$site[i]
  
  pleist <- pleistocene[, sitename]
  holo  <- holocene[, sitename]
  
  # Expand union of species
  pleist_species <- rownames(pleistocene)[pleist > 0]
  holo_species   <- rownames(holocene)[holo > 0]
  union_species  <- union(pleist_species, holo_species)
  
  if (length(union_species) == 0) {
    # no species in either period => turnover undefined or zero
    turnover_df$turnover[i] <- 0
    next
  }
  
  # Intersection, unique to each period
  intersection_species <- intersect(pleist_species, holo_species)
  pleist_only          <- setdiff(pleist_species, holo_species)
  holo_only            <- setdiff(holo_species,   pleist_species)
  
  # Turnover measure
  lost_plus_gained <- length(pleist_only) + length(holo_only)
  turnover_val <- lost_plus_gained / length(union_species)
  
  turnover_df$turnover[i] <- turnover_val
}

#  Merge turnover data with coordinates
coords_turnover <- merge(coords, turnover_df, by = "site")

# Plot a map of the sites colored by turnover
world_map <- subset(map_data("world"), long <= 180)

p_turnover <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "lightgrey", color = NA) +
  geom_point(data = coords_turnover, 
             aes(x = lon, y = lat, color = turnover), 
             size = 2) +
  scale_color_gradient(low = "blue", high = "red", name = "Turnover") +
  coord_map("moll") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 10),
    plot.subtitle = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  ggtitle("Site-specific turnover (Pleistocene→Holocene)")

summary(coords_turnover$turnover)
print(p_turnover)
ggsave("paired_sites_turnover_map.png", p_turnover, width = 6, height = 4, units = "in", dpi = 300, scale = 1)

