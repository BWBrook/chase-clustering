### Quaternary community data ##########################################################################################
{rm(list=ls()); options(scipen=999,digits=9) # clean up and set options
 setwd("C:/Users/bwbrook/Downloads/2_Biogeography_MS")} # Set the working directory to the base directory
#######################################################################################################################

### Import specific package functions #################################################################################
{import::from("dplyr", select, slice)
 import::from("ggplot2", aes, coord_map, geom_point, geom_polygon, ggplot, labs, map_data, theme_void)
 import::from("magrittr", "%>%") }
#######################################################################################################################

### Create custom functions ###########################################################################################
# Compute the fit metric based on the cluster matrix
{computeFitMetric <- function(clusterMatrix) {
  
  clusterSums <- colSums(clusterMatrix) # Sum the columns to get the sum of each cluster
  
  if(sum(clusterSums == 0) > 0) {
    clusterMatrix = clusterMatrix[, clusterSums > 0]
    if(!is.matrix(clusterMatrix)) {
      return(1e6)
    }
    clusterSums = colSums(clusterMatrix)
  } # Remove clusters with zero sum
  
  fitMetric = 0 # Initialize fit metric variable
  
  # Normalize the matrix by dividing each column by its sum
  clusterMatrix = clusterMatrix / clusterSums
  
  # Calculate the sum of each row (i.e., each species across clusters)
  speciesSums = rowSums(clusterMatrix)
  
  # Compute the fit metric
  for (i in 2:ncol(clusterMatrix)) {
    fitMetric = fitMetric + sum((clusterMatrix[, i] * (speciesSums - clusterMatrix[, i]))^0.5)
  }
  
  return(fitMetric)
}

# Create the initial chase matrix
 createChaseMatrix <- function(dataMatrix) {
  # Initialise an empty matrix of appropriate dimensions
  chaseMatrix <- matrix(0, ncol(dataMatrix), ncol(dataMatrix))
  
  # Loop through the columns to populate the chase matrix
  for (i in 2:ncol(dataMatrix)) {
    for (j in 1:(i - 1)) {
      chaseMatrix[i, j] <- sum(dataMatrix[, i] * dataMatrix[, j])
      chaseMatrix[j, i] <- chaseMatrix[i, j]
    }
  }
  
  return(chaseMatrix)
} }
#######################################################################################################################

# Read in the primary datasets
{pleistocene <- read.csv("data/pleistocene_all.csv", row.names = "species") # Pleistocene species-site occ (all)
 holocene <- read.csv("data/holocene_all.csv", row.names = "species") # Holocene sp-site + domesticates
 site_coords <- read.csv("data/all_site_coords.csv", row.names = "site")} # site descriptive variables and attributes

### set the conditions of the analysis ################################################################################
{dataMatrix = pleistocene # data set to test (pleistocene or holocene)
 min_sp = 5 # minimum number of species in a site to include it in the analysis
 max_clusters = 7 # maximum number of clusters to resolve through the chase algorithm
 max_starts = 100 # maximum number of starting configurations, impacts duration of optimisation
 max_shuffles = 10000 } # maximum number of stochastic shuffles to sample when chasing, more has greater stability
#######################################################################################################################

### Set the initiation conditions before the chase clustering algorithm runs ##########################################
{row.names(site_coords) <- as.integer(gsub("s_", "", row.names(site_coords))) # integer site names

 dataMatrix <- dataMatrix %>%
   slice(which(rowSums(.) > 0)) %>% # Remove rows where the sum across the row is <= 0
   select(which(colSums(.) >= min_sp)) %>% # Remove columns where sum < min_sp
   slice(which(rowSums(.) > 0)) # Final filter

 sampleNumber <- colnames(dataMatrix) %>% 
   gsub("s_", "", .) %>% # remove prefix
   as.integer() %>% 
   `[`(., which(colSums(dataMatrix) >= min_sp)) # keep only those meeting min species criteria

 weightVector = colSums(t(t(dataMatrix) / rowSums(dataMatrix))) # used to pick  seed samples in each starting config
 columnWeight = colSums(dataMatrix)  # colSums are used to weight the sample selection step below

 dataMatrix = dataMatrix / columnWeight # columns are standardised by species counts
 chaseMatrix = createChaseMatrix(dataMatrix) # create the chase matrix

 # Initialise the best fit stat and the assignment vector
 finalFitMetric = array(dim = max_clusters, data = 1e6)
 finalClusterAssign = matrix(0, ncol(dataMatrix), max_clusters, dimnames=list(sampleNumber, 1:max_clusters)) }

### Run the chase clustering algorithm runs ###########################################################################
for (numClusters in 2:max_clusters) { # Loop over a range of cluster counts
  
  # Initialize a counter to keep track of the number of improvements found
  improvementsFound = 0
  
  # Loop over a number of different starting configurations
  for (startConfig in 1:max_starts) {
    
    # Initialize cluster assignment vector
    clusterAssign = rep(0, ncol(dataMatrix))
    
    # Randomly select distinct samples to seed the clusters
    clusterAssign[sample(ncol(dataMatrix), numClusters, prob = weightVector)] = 1:numClusters
    
    # Assign the remaining samples to clusters
    for (sampleIdx in which(clusterAssign == 0)) {
      validClusters = which(clusterAssign > 0)
      if (sum(chaseMatrix[sampleIdx, validClusters]) == 0) {
        clusterAssign[sampleIdx] = sample(numClusters, 1)
        next
      }
      clusterAssign[sampleIdx] = clusterAssign[sample(validClusters, 1, prob = chaseMatrix[sampleIdx, validClusters])]
    }
    
    # Initialize the cluster matrix to store the sum of each cluster
    clusterSumMatrix = matrix(0, nrow(dataMatrix), numClusters)
    
    # Populate the cluster sum matrix based on the current cluster assignment
    for (i in 1:ncol(dataMatrix)) {
      clusterSumMatrix[, clusterAssign[i]] = clusterSumMatrix[, clusterAssign[i]] + dataMatrix[, i]
    }
    
    # Initialize the best fit metric for this start configuration
    bestFitMetric = 0
    
    # Loop over a number of shuffles to improve the configuration
    for (shuffleIdx in 1:max_shuffles) {
      
      # Randomly select a sample to move
      selectedSample = sample(1:ncol(dataMatrix), 1)
      oldCluster = clusterAssign[selectedSample]
      
      # Skip if the selected sample would empty its current cluster
      if (sum(clusterAssign == oldCluster) <= 2) {
        next
      }
      
      # Assign the selected sample to a new cluster
      if (sum(chaseMatrix[selectedSample, ]) == 0) {
        newCluster = sample(numClusters, 1)
      } else {
        newCluster = clusterAssign[sample(1:ncol(dataMatrix), 1, prob = chaseMatrix[selectedSample, ])]
      }
      
      # Skip if the sample stays in its current cluster
      if (newCluster == oldCluster) {
        next
      }
      
      # Update the cluster sum matrix to reflect the moved sample
      clusterSumMatrix[, oldCluster] = clusterSumMatrix[, oldCluster] - dataMatrix[, selectedSample]
      clusterSumMatrix[clusterSumMatrix[, oldCluster] < 0, oldCluster] = 0
      clusterSumMatrix[, newCluster] = clusterSumMatrix[, newCluster] + dataMatrix[, selectedSample]
      
      # Compute the new fit metric
      currentFitMetric = computeFitMetric(clusterSumMatrix)
      
      # Update the best fit metric and cluster assignment if improvement is found
      if (currentFitMetric < bestFitMetric || bestFitMetric == 0) {
        bestFitMetric = currentFitMetric
        clusterAssign[selectedSample] = newCluster
        bestClusterAssign = clusterAssign
        bestClusterSumMatrix = clusterSumMatrix
      } else {
        # Revert the changes if no improvement is found
        clusterSumMatrix[, oldCluster] = clusterSumMatrix[, oldCluster] + dataMatrix[, selectedSample]
        clusterSumMatrix[, newCluster] = clusterSumMatrix[, newCluster] - dataMatrix[, selectedSample]
        clusterSumMatrix[clusterSumMatrix[, newCluster] < 0, newCluster] = 0
      }
    }
    
    # Save the best fit metric and cluster assignment for this number of clusters
    if (bestFitMetric < finalFitMetric[numClusters] - 1e-8) {
      finalFitMetric[numClusters] = bestFitMetric
      finalClusterAssign[, numClusters] = bestClusterAssign
      improvementsFound = 1
      
      cat('\rfound new ', numClusters, startConfig, bestFitMetric)
      plot(site_coords[as.character(sampleNumber),1], site_coords[as.character(sampleNumber),2], 
           col=hsv(h=(bestClusterAssign-1)^0.5/numClusters^0.5), pch=bestClusterAssign,
           xlab="Lon", ylab="Lat", main="Current best geographic chase clustering map")

      # If the fit metric is nearly the same as the best found so far
    } else if (bestFitMetric > finalFitMetric[numClusters] - 1e-8 && bestFitMetric < finalFitMetric[numClusters] + 1e-8) {
      # Increment the counter for the number of times this fit metric has been found
      improvementsFound = improvementsFound + 1
      cat('\rfound again    ', improvementsFound)
    }
    
    # Display the progress
    cat('\n', startConfig, bestFitMetric)
  }
  
  # Display the best fit metric and the number of times it was found for this number of clusters
  cat('\nbest', numClusters, improvementsFound, finalFitMetric[numClusters], '\n')
}
#######################################################################################################################
