### Quaternary community data

setwd("C:/git/chase-clustering") # choose appropriate location

## Import specific package functions
import::from("src/qm_clustchase_func.r", .all=T) # import custom functions

# --- Main Script ---
# Read primary datasets
pleistocene <- read.csv("data/pleistocene_all.csv", row.names = "species")   # Pleistocene species-site occupancy
holocene    <- read.csv("data/holocene_wild.csv", row.names = "species")      # Holocene wild species-site (without domesticates)
#holocene    <- read.csv("data/holocene_all.csv", row.names = "species")      # Holocene all species-site (with domesticates)
site_coords <- read.csv("data/all_site_coords.csv", row.names = "site")       # Site descriptive variables and coordinates

# Set analysis parameters
dataMatrix   <- holocene # Choose dataset: pleistocene or holocene
alpha        <- -0.1        # Set to negative (e.g., -0.1) to encourage more clustering, 0 for no penalization, positive (e.g., 0.1) for fewer clusters
min_sp       <- 5           # Minimum species count per site to include
max_clusters <- 15          # Maximum number of clusters to test
max_starts   <- 250         # Number of starting configurations (affects run duration and max clusters found)
max_shuffles <- 10000       # Number of stochastic shuffles for optimization (affects stability of stochastic search)

preproc <- prepareChaseClusteringData(dataMatrix = dataMatrix, min_sp = min_sp, max_clusters = max_clusters)

### Run the chase clustering algorithm runs
result <- runChaseClustering(preproc, max_clusters = max_clusters, max_starts = max_starts, max_shuffles = max_shuffles, 
                             site_coords = site_coords, alpha = alpha, verbose = TRUE)

print(result$finalFitMetric); print(result$universalBestClusters)

dat_clust <- result$universalBestAssign %>%
  enframe(name = "site", value = "cluster") %>%  # Convert named vector to a tibble
  mutate(site = paste0("s_", site)) %>%          # Prepend "s_" to match site_coords row names
  left_join(rownames_to_column(site_coords, var = "site"), by = "site")

# Plot clustering results (dat_clust must have lon, lat, and cluster columns)
p_chase <- plot_world_map(point_data = dat_clust,
                         color_var = "cluster",
                         point_size = 2,
                         legend_position = "none",
                         title = NULL, 
                         subtitle = NULL,
                         save_file = 'Holocene_all_chase_a-01.png')
print(p_chase)
