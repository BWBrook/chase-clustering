### Quaternary community data ##########################################################################################
{rm(list=ls()); options(scipen=999,digits=9) # clean up and set options
 setwd("C:/Users/bwbrook/Downloads/2_Biogeography_MS")} # Set the working directory to the base directory
#######################################################################################################################

### Import specific package functions #################################################################################
{import::from("ClustGeo", choicealpha, hclustgeo)
 import::from("dplyr", across, filter, mutate, select, ungroup)
 import::from("ggplot2", aes, coord_map, geom_point, geom_polygon, ggplot, labs, map_data, theme_void)
 import::from("magrittr", "%>%")
 import::from("sf", st_as_sf, st_distance)}
#######################################################################################################################

# Read in the primary datasets
{pleistocene <- read.csv("data/pleistocene_all.csv", row.names = "species") # Pleistocene species-site occ (all)
 holocene <- read.csv("data/holocene_all.csv", row.names = "species") # Holocene sp-site + domesticates
 site_coords <- read.csv("data/all_site_coords.csv", row.names = "site")} # site descriptive variables and attributes

### set the conditions of the analysis ################################################################################
{dat = pleistocene # data set to test (pleistocene or holocene)
 min_sp = 5 # minimum number of species in a site to include it in the analysis
 k = 7} # set number of clusters
#######################################################################################################################

{dat <- dat %>% # Create filtered occupancy or count dataframe for sites (species composition)
    select(-species) %>%
    select(where(~sum(.x != 0) >= min_sp)) %>%
    mutate(across(everything(), ~. / sum(.))) # create proportional matrix

 # Dataframe to store the site, coordinates and cluster identity
 dat_clust <- site_coords %>%
   filter(site %in% colnames(dat)) %>%
   select(site, lon, lat) %>%
   mutate(cluster = NA_real_) }

{d0 <- dist(t(dat)) # create a distance matrix based on the occupancy data

 # calculate the great-circle distance matrix on the geographic coordinates
 d1 <- dat_clust %>%
   st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # WGS84
   st_distance() %>%
   as.dist() }

# Ward Hierarchical clustering based on distance matrix of features mixed with geographic coordinates
{tree <- hclustgeo(d0)
 dat_clust$cluster <- cutree(tree, k)} # Cut the tree at k to set clusters

# Examine plot of dendrogram with selected k clusters
{plot(tree, hang = -1, label = FALSE, xlab = "", sub = "", main = "Ward dendrogram with feature only")
 rect.hclust(tree, k = k, border = c(1:k))
 legend("topright", legend = paste("cluster", 1:k), fill = 1:k, bty = "n", border = "white")}

# use the choicealpha() function to choose an optimal alpha value (0 = pure site, 1 = pure geog)
{cr <- choicealpha(d0, d1, seq(0, 1, by = 0.05), k, graph=T)
 alpha <- cr$range.alpha[nrow(cr$Qnorm) - which.max(rev(rowSums(cr$Qnorm))) + 1]} # max inertia var expl
#alpha<-0}

# perform spatially constrained hierarchical clustering, with intertia-optimised alpha mixing parameter
{dat_clust$cluster <- cutree(hclustgeo(d0, d1, alpha), k)
 
 world_map = map_data("world") %>% 
    filter(! long > 180)

 ggplot() + # Use ggplot to create the map with the different clusters coloured
   geom_polygon(data=world_map, aes(x=long, y=lat, group=group), fill="lightgrey", color=NA) +
   geom_point(data=dat_clust, aes(x=lon, y=lat, color=factor(cluster)), size=4) +
   coord_map("moll") +
   theme_void() +
   labs(color='Cluster') }
#######################################################################################################################
