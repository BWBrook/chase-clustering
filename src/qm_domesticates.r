# Identifies domesticated species in the Holocene data to produce:
# (1) A summary of domesticates' frequency across sites
# (2) An example "chase clustering" for holocene_all vs. holocene_wild to see how domesticates affect cluster assignments
# (3) A table of results highlighting the impact of domesticates

setwd("C:/git/chase-clustering")  # adapt as needed

import::from("dplyr", summarise, filter)
import::from("fossil", adj.rand.index)
import::from("src/qm_clustchase_func.r", .all=T) # import custom functions

# Read primary datasets
holocene_all    <- read.csv("data/holocene_all.csv",    row.names = 1, check.names = FALSE)
holocene_wild   <- read.csv("data/holocene_wild.csv",   row.names = 1, check.names = FALSE)

# Identify domesticated species by row names present in holocene_all but absent in holocene_wild
domesticated_species <- setdiff(rownames(holocene_all), rownames(holocene_wild))

# Summaries
cat("Number of total Holocene species (all):", nrow(holocene_all), "\n")
cat("Number of Holocene wild species:", nrow(holocene_wild), "\n")
cat("Number of Holocene domesticated species:", length(domesticated_species), "\n")

# Summarise how many sites (columns) contain at least one domesticated species
holocene_all_mat <- as.matrix(holocene_all)  # ensure numeric or logical
dom_matrix <- holocene_all_mat[domesticated_species, , drop = FALSE]
sites_with_dom <- apply(dom_matrix, 2, function(x) any(x > 0))
num_sites_with_dom <- sum(sites_with_dom)
cat("Number of sites with >= 1 domesticated species:", num_sites_with_dom, "\n")

# Compute mean # of domesticated species per site
mean_dom_per_site <- mean(colSums(dom_matrix > 0))
cat("Mean # of domesticated species per site (across all Holocene sites):", 
    round(mean_dom_per_site, 2), "\n")

# Set analysis parameters
alpha        <- -0.1        # Set to negative (e.g., -0.1) to encourage more clustering, 0 for no penalization, positive (e.g., 0.1) for fewer clusters
min_sp       <- 5           # Minimum species count per site to include
max_clusters <- 15          # Maximum number of clusters to test
max_starts   <- 250         # Number of starting configurations (affects run duration and max clusters found)
max_shuffles <- 10000       # Number of stochastic shuffles for optimization (affects stability of stochastic search)

# (a) Holocene All
prep_all <- prepareChaseClusteringData(dataMatrix = holocene_all, 
                                       min_sp = min_sp, 
                                       max_clusters = max_clusters)

res_all  <- runChaseClustering(prep_all, 
                               max_clusters = max_clusters, 
                               max_starts = max_starts, 
                               max_shuffles = max_shuffles,
                               site_coords = NULL,
                               alpha = alpha,
                               verbose = TRUE)

k_best_all <- which.min(res_all$finalFitMetric[2:max_clusters]) + 1
clusters_all <- res_all$finalClusterAssign[, k_best_all]

# (b) Holocene Wild Only
prep_wild <- prepareChaseClusteringData(dataMatrix = holocene_wild,
                                        min_sp = min_sp,
                                        max_clusters = max_clusters)
res_wild  <- runChaseClustering(prep_wild,
                                max_clusters = max_clusters,
                                max_starts = max_starts,
                                max_shuffles = max_shuffles,
                                site_coords = NULL,
                                alpha = alpha,
                                verbose = TRUE)

k_best_wild <- which.min(res_wild$finalFitMetric[2:max_clusters]) + 1
clusters_wild <- res_wild$finalClusterAssign[, k_best_wild]

# Restrict to the shared site set
shared_sites <- intersect(names(clusters_all), names(clusters_wild))
common_clust_all  <- clusters_all[shared_sites]
common_clust_wild <- clusters_wild[shared_sites]

# Compute how many sites changed cluster membership
num_changed_sites <- sum(common_clust_all != common_clust_wild)
cat("Number of sites (in common) that differ in cluster membership after removing domesticates:", 
    num_changed_sites, "\n")

# Compute an Adjusted Rand Index for the partition difference
ari_value <- fossil::adj.rand.index(common_clust_all, common_clust_wild)
cat("Adjusted Rand Index (ARI) between Holocene All vs. Wild-only cluster partitions:", 
    round(ari_value, 3), "\n")

# Create a Summary Table
res_table <- data.frame(
  TotalHoloceneSpecies = nrow(holocene_all),
  DomesticatedSpecies  = length(domesticated_species),
  SitesWithDomestics   = num_sites_with_dom,
  MeanDomesticsPerSite = round(mean_dom_per_site, 2),
  BestK_all            = k_best_all,
  BestK_wild           = k_best_wild,
  ChangedSites         = num_changed_sites,
  ARI                  = round(ari_value, 3)
)

print(res_table)
write.csv(res_table, "domestics_impact_summary.csv", row.names = FALSE)

## For plot in Supplementary Information

## Import specific package functions
import::from("src/qm_clustchase_func.r", .all=T) # import custom functions
import::from("dplyr", filter)
library(ggplot2)
library(sf)

# Read Holocene data
holocene_all <- read.csv("data/holocene_all.csv", row.names = 1, check.names = FALSE)
holocene_wild <- read.csv("data/holocene_wild.csv", row.names = 1, check.names = FALSE)

# Read site coordinates (ensure this CSV has columns: "site", "lon", "lat")
site_coords <- read.csv("data/all_site_coords.csv")

# Identify domesticated species
holocene_all_mat <- as.matrix(holocene_all)
domesticated_species <- setdiff(rownames(holocene_all), rownames(holocene_wild))
dom_matrix <- holocene_all_mat[domesticated_species, , drop = FALSE]

# Calculate proportion of domesticated species per site
total_species_per_site <- colSums(holocene_all_mat > 0)
domestics_per_site <- colSums(dom_matrix > 0)
prop_domestics_per_site <- domestics_per_site / total_species_per_site

# Merge with coordinates
prop_dom_df <- data.frame(site = names(prop_domestics_per_site), proportion_domestics = prop_domestics_per_site)
coords_domestics <- merge(site_coords, prop_dom_df, by = "site")

# Plot map of proportion domesticated
world_map <- subset(map_data("world"), long <= 180)

p_domestics <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "lightgrey", color = NA) +
  geom_point(data = coords_domestics, 
             aes(x = lon, y = lat, color = proportion_domestics), 
             size = 2) +
  scale_color_gradient(low = "blue", high = "red", name = "Proportion domesticated") +
  coord_map("moll") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 10),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  ggtitle("Proportion of Domesticated Species per Holocene Site")

summary(coords_domestics$proportion_domestics)
print(p_domestics)
ggsave("supplement_proportion_domestics_map.png", p_domestics, width = 6, height = 4, units = "in", dpi = 300, scale = 1)
