# Hierarchical Geographic Clustering Using Ward's Method with Spatial Constraints
# Quaternary Community Analysis

# Load necessary functions from packages
import::from("ClustGeo", choicealpha, hclustgeo)
import::from("dplyr", across, filter, mutate, select, ungroup)
import::from("magrittr", "%>%")
import::from("sf", st_as_sf, st_distance)
import::from("tibble", rownames_to_column)
library(ggplot2)

setwd("C:/git/chase-clustering") # choose appropriate location

# Read in primary datasets (ensure the working directory is set appropriately)
pleistocene <- read.csv("data/pleistocene_all.csv", row.names = "species")  # Pleistocene species-site occupancy
holocene    <- read.csv("data/holocene_all.csv", row.names = "species")     # Holocene species-site (with domesticates)
site_coords <- read.csv("data/all_site_coords.csv", row.names = "site")      # Site descriptive variables and coordinates

# Set analysis parameters
dat    <- pleistocene       # Data set to test (pleistocene or holocene)
min_sp <- 5                 # Minimum number of species per site
k      <- 7                 # Desired number of clusters

# Create filtered occupancy matrix:
# - Exclude sites with fewer than min_sp species.
# - Convert counts to proportions (each site's column sums to 1).
dat <- dat %>%
  select(where(~ sum(.x != 0) >= min_sp)) %>%# Filter columns by species count
  mutate(across(everything(), ~ . / sum(.))) # Standardize each column

# Prepare a site coordinate dataframe and initialize the cluster column
dat_clust <- site_coords %>%
  rownames_to_column(var = "site") %>%
  filter(site %in% colnames(dat)) %>%
  select(site, lon, lat) %>%
  mutate(cluster = NA_real_)

# Compute dissimilarity matrix (feature space) based on species occupancy data
d0 <- dist(t(dat))

# Compute geographic distance matrix (constraint space)
d1 <- dat_clust %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%  # Ensure coordinates are in WGS84
  st_distance() %>%                                   # Compute pairwise distances (m)
  as.dist()

# Initial hierarchical clustering using Ward's method on features only
tree <- hclustgeo(d0)
dat_clust$cluster <- cutree(tree, k)  # Cut tree to assign clusters

# Plot dendrogram for feature-only clustering
plot(tree, hang = -1, label = FALSE, xlab = "", sub = "", main = "Ward Dendrogram (Feature Only)")
rect.hclust(tree, k = k, border = 1:k)
legend("topright", legend = paste("Cluster", 1:k), fill = 1:k, bty = "n", border = "white")

# Determine optimal alpha by maximizing variance explanation:
alpha_seq <- seq(0, 1, by = 0.01)  # fine-scale sampling
cr <- choicealpha(d0, d1, alpha_seq, k, graph = TRUE)
df_qnorm <- data.frame(alpha = cr$range.alpha, sumQ = rowSums(cr$Qnorm))
candidate_alphas <- df_qnorm$alpha[abs(df_qnorm$sumQ - max(df_qnorm$sumQ)) < 1e-6]
alpha_opt <- max(candidate_alphas); alpha_opt

# Perform spatially constrained hierarchical clustering using the optimal alpha
tree_spatial <- hclustgeo(d0, d1, alpha_opt)
dat_clust$cluster <- cutree(tree_spatial, k)

# Plot the clustering results on a world map
world_map <- subset(map_data("world"), long <= 180)

p_ward <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgrey", color = NA) +
  geom_point(data = dat_clust, aes(x = lon, y = lat, color = factor(cluster)), size = 2) +
  coord_map("moll") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Display and save individual plot to PNG files
p_ward # Note: designed to be included in panel plot with additional title and labels
ggsave("Pleistocene_Ward.png", p_ward, width = 6, height = 4, units = "in", dpi = 300, scale = 1)

