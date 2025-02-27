# qe_clustering and plotting functions

import::here("dplyr", left_join, mutate, select, slice)
import::here("magrittr", "%>%")
import::here("tibble", enframe, rownames_to_column)
import::here("ggplot2", ggplot, geom_polygon, geom_point, aes, aes_string, coord_map, theme_void, theme, element_text, labs, ggsave, map_data, margin)

# Preprocessing: ensure numeric, remove zero-row species, filter columns by min_sp, build chaseMatrix
prepareChaseClusteringData <- function(dataMatrix, min_sp, max_clusters) {
  # dataMatrix is read from CSV, so it may have come in as data.frame with rownames as species
  # Make sure everything except rownames is numeric
  dataMatrix <- as.matrix(dataMatrix)  # tries to convert to numeric matrix; errors if any non-numeric cell
  
  # Remove rows that are entirely zero (no presence across all sites)
  rowTotals <- rowSums(dataMatrix)
  dataMatrix <- dataMatrix[rowTotals > 0, , drop = FALSE]
  
  # Filter columns (sites) to only those with colSum >= min_sp
  colTotals <- colSums(dataMatrix)
  keepCols  <- which(colTotals >= min_sp)
  dataMatrix <- dataMatrix[, keepCols, drop = FALSE]
  
  # Possibly remove species rows that became zero after col filtering
  rowTotals <- rowSums(dataMatrix)
  dataMatrix <- dataMatrix[rowTotals > 0, , drop = FALSE]
  
  # sampleNumber is derived from the column names if they are like "s_1", "s_4", etc.
  # Suppose we want them as integers. For safety, we can do:
  colNames <- colnames(dataMatrix)
  sampleNumber <- gsub("^s_", "", colNames)   # strip off "s_"
  sampleNumber <- as.integer(sampleNumber)    # convert to integer
  
  # For weighting the random seeds, you could use column sums or just uniform:
  weightVector <- colSums(dataMatrix)
  
  # Now build the chaseMatrix via crossprod: size is (#columns) x (#columns)
  chaseMatrix <- crossprod(dataMatrix)
  
  # Initialize final containers
  finalFitMetric <- rep(1e6, max_clusters)
  finalClusterAssign <- matrix(0, ncol(dataMatrix), max_clusters,
                               dimnames = list(sampleNumber, paste0("C", 1:max_clusters)))
  
  return(list(
    dataMatrix         = dataMatrix,
    sampleNumber       = sampleNumber,
    weightVector       = weightVector,
    chaseMatrix        = chaseMatrix,
    finalFitMetric     = finalFitMetric,
    finalClusterAssign = finalClusterAssign
  ))
}

# A simple pairwise overlap measure: bigger overlap => bigger metric => "worse" => we want to minimize
computeFitMetric <- function(clusterMatrix) {
  clusterSums <- colSums(clusterMatrix)
  
  # Remove clusters with zero total
  if (any(clusterSums == 0)) {
    clusterMatrix <- clusterMatrix[, clusterSums > 0, drop = FALSE]
    # If that leaves < 2 clusters, return a large penalty
    if (ncol(clusterMatrix) < 2) {
      return(1e6)
    }
    clusterSums <- colSums(clusterMatrix)
  }
  
  # Normalize columns
  clusterMatrix <- sweep(clusterMatrix, 2, clusterSums, FUN = "/")
  
  # Pairwise overlap
  nclust <- ncol(clusterMatrix)
  overlap <- 0
  for (i in seq_len(nclust - 1)) {
    for (j in (i + 1):nclust) {
      overlap_ij <- sum(sqrt(clusterMatrix[, i] * clusterMatrix[, j]))
      overlap <- overlap + overlap_ij
    }
  }
  return(overlap)
}

# Main chase clustering. We minimize the (overlap + alpha * numClusters). See Fig. 1 flowchart in main paper.
runChaseClustering <- function(preproc,
                               max_clusters,
                               max_starts,
                               max_shuffles,
                               site_coords,
                               alpha   = 0,     # set > 0 to penalize more clusters
                               verbose = TRUE) {
  
  dataMatrix         <- preproc$dataMatrix
  sampleNumber       <- preproc$sampleNumber
  weightVector       <- preproc$weightVector
  chaseMatrix        <- preproc$chaseMatrix
  finalFitMetric     <- preproc$finalFitMetric
  finalClusterAssign <- preproc$finalClusterAssign
  
  universalBestFit      <- 1e6
  universalBestClusters <- NA
  universalBestAssign   <- NULL
  
  noImprovementCount <- 0
  
  for (numClusters in 2:max_clusters) {
    improvementsFound <- 0
    
    for (startConfig in seq_len(max_starts)) {
      # Initialize cluster assignment
      clusterAssign <- integer(ncol(dataMatrix))
      
      # Randomly seed the clusters
      seedIndices <- sample(ncol(dataMatrix), numClusters, prob = weightVector)
      clusterAssign[seedIndices] <- seq_len(numClusters)
      
      # Assign remaining samples
      for (sampleIdx in which(clusterAssign == 0)) {
        validClusters <- which(clusterAssign > 0)
        rowProbs <- chaseMatrix[sampleIdx, validClusters]
        if (sum(rowProbs) == 0) {
          clusterAssign[sampleIdx] <- sample.int(numClusters, 1)
        } else {
          chosenCol <- sample(validClusters, 1, prob = rowProbs)
          clusterAssign[sampleIdx] <- clusterAssign[chosenCol]
        }
      }
      
      # Build cluster sum matrix (rows=species, cols=clusters)
      clusterSumMatrix <- matrix(0, nrow(dataMatrix), numClusters)
      for (colIdx in seq_len(ncol(dataMatrix))) {
        clusterSumMatrix[, clusterAssign[colIdx]] <- clusterSumMatrix[, clusterAssign[colIdx]] + dataMatrix[, colIdx]
      }
      
      bestFitMetric <- Inf
      bestClusterAssign <- clusterAssign
      bestClusterSumMatrix <- clusterSumMatrix
      
      # Shuffle to improve
      for (shuffleIdx in seq_len(max_shuffles)) {
        selectedSample <- sample.int(ncol(dataMatrix), 1)
        oldCluster <- clusterAssign[selectedSample]
        
        # Avoid emptying that cluster if only 1 sample is in it
        if (sum(clusterAssign == oldCluster) <= 1) next
        
        rowProbs <- chaseMatrix[selectedSample, ]
        if (sum(rowProbs) == 0) {
          newCluster <- sample.int(numClusters, 1)
        } else {
          chosenCol  <- sample.int(ncol(dataMatrix), 1, prob = rowProbs)
          newCluster <- clusterAssign[chosenCol]
        }
        
        # If no change, skip
        if (newCluster == oldCluster) next
        
        # Move sample
        clusterSumMatrix[, oldCluster] <- clusterSumMatrix[, oldCluster] - dataMatrix[, selectedSample]
        clusterSumMatrix[, newCluster] <- clusterSumMatrix[, newCluster] + dataMatrix[, selectedSample]
        
        # Compute penalized fit
        rawFit      <- computeFitMetric(clusterSumMatrix)
        if(rawFit==0) rawFit <- 1e6 # skip degenerate cases where rawFit == 0 (typically with clusters == 2)
        penalizedFit <- rawFit + alpha * numClusters
        
        if (penalizedFit < bestFitMetric) {
          bestFitMetric <- penalizedFit
          clusterAssign[selectedSample] <- newCluster
          bestClusterAssign <- clusterAssign
          bestClusterSumMatrix <- clusterSumMatrix
        } else {
          # Revert
          clusterSumMatrix[, newCluster] <- clusterSumMatrix[, newCluster] - dataMatrix[, selectedSample]
          clusterSumMatrix[, oldCluster] <- clusterSumMatrix[, oldCluster] + dataMatrix[, selectedSample]
        }
      }
      
      # Compare to finalFitMetric
      if (bestFitMetric < finalFitMetric[numClusters] - 1e-8) {
        finalFitMetric[numClusters] <- bestFitMetric
        finalClusterAssign[, numClusters] <- bestClusterAssign
        improvementsFound <- improvementsFound + 1
        
        if (verbose) {
          cat(sprintf("\nFound new best for cluster=%d start=%d fit=%.6f\n", 
                      numClusters, startConfig, bestFitMetric))
        }
        # Optional plot
        #idx <- paste0("s_", sampleNumber)
        #plot_coords <- site_coords[idx, , drop = FALSE]
        #lon <- as.numeric(as.character(plot_coords[, 1]))
        #lat <- as.numeric(as.character(plot_coords[, 2]))
        #if (length(lon) && length(lat) && all(is.finite(lon)) && all(is.finite(lat))) {
        #  plot(lon, lat,
        #       col = hsv(h = (bestClusterAssign - 1)^0.5 / numClusters^0.5),
        #       pch = bestClusterAssign,
        #       xlab = "Lon", ylab = "Lat",
        #       main = sprintf("Best for %d clusters (current)", numClusters))
        #}
      } else if (abs(bestFitMetric - finalFitMetric[numClusters]) < 1e-8) {
        improvementsFound <- improvementsFound + 1
      }
      
      if (verbose) {
        cat(sprintf("StartConfig %d final penalized fit=%.6f\n", startConfig, bestFitMetric))
      }
    }
    
    if (verbose) {
      cat(sprintf("NumClusters=%d improved=%d bestFit=%.6f\n", 
                  numClusters, improvementsFound, finalFitMetric[numClusters]))
    }
    
    # Update universal best
    if (finalFitMetric[numClusters] < universalBestFit) {
      universalBestFit      <- finalFitMetric[numClusters]
      universalBestClusters <- numClusters
      universalBestAssign   <- finalClusterAssign[, numClusters]
      
      # Optional plot
      idx <- paste0("s_", sampleNumber)
      plot_coords <- site_coords[idx, , drop = FALSE]
      lon <- as.numeric(as.character(plot_coords[, 1]))
      lat <- as.numeric(as.character(plot_coords[, 2]))
      if (length(lon) && length(lat) && all(is.finite(lon)) && all(is.finite(lat))) {
        plot(lon, lat,
             col = hsv(h = (universalBestAssign - 1)^0.5 / universalBestClusters^0.5),
             pch = universalBestAssign,
             xlab = "Lon", ylab = "Lat",
             main = sprintf("Universal Best: %d clusters, penalizedFit=%.6f", 
                            universalBestClusters, universalBestFit))
      }
      noImprovementCount <- 0
    } else {
      noImprovementCount <- noImprovementCount + 1
      if (noImprovementCount > 2) {
        if (verbose) cat("No improvement for +3 cluster sizes, stopping.\n")
        break
      }
    }
  }
  
  return(list(
    finalFitMetric       = finalFitMetric,
    finalClusterAssign   = finalClusterAssign,
    universalBestFit     = universalBestFit,
    universalBestClusters= universalBestClusters,
    universalBestAssign  = universalBestAssign
  ))
}

plot_world_map <- function(point_data, 
                           color_var = "cluster", 
                           point_size = 2, 
                           legend_position = "none", 
                           title = NULL, 
                           subtitle = NULL, 
                           plot_margin = margin(10, 10, 10, 10), 
                           save_file = NULL, 
                           width = 6, height = 4, 
                           units = "in", dpi = 300, scale = 1) {

  p <- ggplot() +
    geom_polygon(data = subset(map_data("world"), long <= 180), 
                 aes(x = long, y = lat, group = group), 
                 fill = "lightgrey", color = NA) +
    geom_point(data = point_data, 
               aes_string(x = "lon", y = "lat", 
                          color = paste0("factor(", color_var, ")")), 
               size = point_size) +
    coord_map("moll") +
    theme_void() +
    theme(legend.position = legend_position,
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.margin = plot_margin) +
    labs(title = title, subtitle = subtitle)
  
  if (!is.null(save_file)) {
    ggsave(filename = save_file, plot = p, 
           width = width, height = height, units = units, dpi = dpi, scale = scale)
  }
  
  return(p)
}

