# domesticates_analysis.R

# This script identifies domesticated species in the Holocene data by comparing
# row names of 'holocene_all.csv' to 'holocene_wild.csv'. It then produces:
# (1) A summary of domesticates' frequency across sites
# (2) An example "chase clustering" for holocene_all vs. holocene_wild 
#     to see how domesticates affect cluster assignments
# (3) A table of results highlighting the impact of domesticates

setwd("C:/git/chase-clustering")  # adapt as needed

import::from("dplyr", summarise, mutate, select, filter)
import::from("fossil", adj.rand.index)
import::from("magrittr", "%>%")
import::from("src/qe_biogr_clustchase_func.r", .all=T) # import custom functions
library(ggplot2)

#-----------------------
# 1. Load data
#-----------------------
pleistocene_all <- read.csv("data/pleistocene_all.csv", row.names = 1, check.names = FALSE)
holocene_all    <- read.csv("data/holocene_all.csv",    row.names = 1, check.names = FALSE)
holocene_wild   <- read.csv("data/holocene_wild.csv",   row.names = 1, check.names = FALSE)

# Identify domesticated species by row names present in holocene_all but absent in holocene_wild
domesticated_species <- setdiff(rownames(holocene_all), rownames(holocene_wild))

#-----------------------
# 2. Summaries
#-----------------------
cat("Number of total Holocene species (all):", nrow(holocene_all), "\n")
cat("Number of Holocene wild species:", nrow(holocene_wild), "\n")
cat("Number of Holocene domesticated species:", length(domesticated_species), "\n")

# Summarize how many sites (columns) contain at least one domesticated species
# We treat these CSVs as presence/absence or counts > 0 => presence
holocene_all_mat <- as.matrix(holocene_all)  # ensure numeric or logical
dom_matrix       <- holocene_all_mat[domesticated_species, , drop = FALSE]

# Sites that have at least one domesticate = any row is > 0
sites_with_dom <- apply(dom_matrix, 2, function(x) any(x > 0))
num_sites_with_dom <- sum(sites_with_dom)
cat("Number of sites with >= 1 domesticated species:", num_sites_with_dom, "\n")

# Compute mean # of domesticated species per site
mean_dom_per_site <- mean(colSums(dom_matrix > 0))
cat("Mean # of domesticated species per site (across all Holocene sites):", 
    round(mean_dom_per_site, 2), "\n")

#-----------------------
# 3. Example: Effect on Clustering
#-----------------------
# We'll do a short chase clustering for holocene_all vs. holocene_wild 
# to see if removing domestics changes cluster membership significantly.

alpha        <- -0.2
max_clusters  <- 10
max_starts    <- 10
max_shuffles  <- 1000
min_sp        <- 5  # minimum species per site

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
                               verbose = FALSE)

# We'll pick the "best" cluster assignment for a certain K or the universal best
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
                                verbose = FALSE)

k_best_wild <- which.min(res_wild$finalFitMetric[2:max_clusters]) + 1
clusters_wild <- res_wild$finalClusterAssign[, k_best_wild]

# Since the site sets might differ if some columns got dropped 
# (due to min_sp in wild-only data), let's restrict to the shared site set:
shared_sites <- intersect(names(clusters_all), names(clusters_wild))
common_clust_all  <- clusters_all[shared_sites]
common_clust_wild <- clusters_wild[shared_sites]

# Compute how many sites changed cluster membership
# We'll do a simple count of sites that differ
num_changed_sites <- sum(common_clust_all != common_clust_wild)
cat("Number of sites (in common) that differ in cluster membership after removing domesticates:", 
    num_changed_sites, "\n")

# Compute an Adjusted Rand Index for the partition difference
ari_value <- fossil::adj.rand.index(common_clust_all, common_clust_wild)
cat("Adjusted Rand Index (ARI) between Holocene All vs. Wild-only cluster partitions:", 
    round(ari_value, 3), "\n")

#-----------------------
# 4. Create a Summary Table
#-----------------------
# Suppose we want a small table summarizing:
#   1) total Holocene species (all)
#   2) # domesticated species
#   3) # sites with domestics
#   4) mean # domestics per site
#   5) best K in all
#   6) best K in wild
#   7) # changed sites
#   8) ARI

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

# You can write it to CSV if you like:
write.csv(res_table, "domestics_impact_summary.csv", row.names = FALSE)
