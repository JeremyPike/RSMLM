
library("dbscan")
library("deldir")
library("igraph")

clusterCCs <- function(coords, r) {
  
  # fixed radius neighbour search
  nn <- frNN(coords, eps = r)
  # create igraph object from veronoi triangulation
  g <- graph_from_adj_list(adjacencylist(nn))
  # find the connected components (clusters of the graph)
  ccs <- components(g)
  # extract cluster indices
  clusterIndices <- ccs$membership
  
  # set noise (clusters with one detecion) to zero and reindex remaining clusters
  numClusters <- length(ccs$csize)
  count <- 0
  for (i in 1 : numClusters) {
    if (ccs$csize[i] == 1) {
      clusterIndices[clusterIndices == i] <- 0
    } else {
      count <- count + 1
      clusterIndices[clusterIndices == i] <- count
    }
  }
  
  return(clusterIndices)
  
}

clusterTomato <- function(coords, r, threshold) {
  
  # compute density using number of other detections within fixed search radius
  density <- pointdensity(coords, r, type = "frequency")
  clusterIndices <- tomatoDens(coords, density, r, threshold)$clusters
  return(clusterIndices)
}

clusterDBSCAN <- function(coords, r, minPoints) {
  
  # perform DBSCAN using "dbscan" pacakage
  fr <- frNN(coords, r)
  clusterIndices <- dbscan(fr, minPts = minPoints, borderPoints=TRUE)$cluster
  
  return(clusterIndices)
  
}
clusterVoronoi <- function(coords, threshold) {
  
  if (sum(duplicatedxy(coords[, 1], coords[, 2])) > 0) {
    stop("Detection list should not contain duplicated points")
  } else {
    numDetections <- dim(coords)[1]
    # compute veronoi tesselation
    vtess <- deldir(coords[, 1], coords[, 2])
    tileAreas <- vtess$summary$dir.area
   
    
    # create igraph object from veronoi triangulation
    g <- make_empty_graph(n = numDetections, directed = FALSE)
    g <- add_edges(g, c(rbind(vtess$delsgs$ind1, vtess$delsgs$ind2)))
    
    # compute first rank density
    density <- c()
    neigh <- adjacent_vertices(g, seq(1, numDetections, by = 1))
    for (i in 1 : numDetections) {
      density[i] = (1 + length(neigh[[i]])) / (tileAreas[i] + sum(tileAreas[neigh[[i]]]))
    }
    # find tiles with density larger than threshold
    tileIndices <- density > threshold
    
    # delete tiles with area larger than specified maxiumum
    g <- delete_vertices(g, which(!tileIndices))
    # compute the connected components (clusters of the graph)
    
    ccs <- components(g)
    # extract cluster indices as vector with zero value representing no cluster assignment
    clusterIndices <- rep(0, numDetections)
    clusterIndices[which(tileIndices)] <- ccs$membership
    
    return(clusterIndices)
  }
}

clusterRipley <- function(coords, r, threshold, ROIArea) {
  
  numDetections <- dim(coords)[1]
  
  threeD <- FALSE
  if (dim(coords)[2] == 3) {
    threeD <- TRUE
  }
  fr <- frNN(coords, eps = r)
  
  numNeighbours <- c()
  for (k in 1 : numDetections) {
    # dor each detection calculate number of other detections within search radius
    numNeighbours[k] <- length(fr$id[[k]])
  }
  
  # calculate ripleys K and L - r functions
  K <- numNeighbours * ROIArea / (numDetections - 1)
  
  if (threeD) {
    LminR <- (3 * K / pi / 4)^(1/3) - r
  } else {
    LminR <- sqrt(K / pi) - r
  }
  
  clusterIndices <- rep(0, numDetections)
  filteredDetections <- LminR > threshold
  numFilteredDetections <- sum(filteredDetections)
  
  if(numFilteredDetections == 1) {

    clusterIndices[filteredDetections] <- 1
    
  } else if (numFilteredDetections > 1)  {
    g <- graph_from_adj_list((adjacencylist)(fr))
    g <- delete_vertices(g, which(!filteredDetections)) 
    ccs <- components(g)
    
    # extract cluster indices as vector with zero value representing no cluster assignment
    
    clusterIndices[filteredDetections] <- ccs$membership
  }
  return(clusterIndices)
}

tomatoDiagram <- function(coords, r) {
  
  # compute density using number of other detections within fixed search radius
  density <- pointdensity(coords, r, type = "frequency")
  diagram <- tomatoDens(coords, density, r, 1e20)$diagram
  return(diagram)
}

filtClustDetections <- function(clusterIndices, minDetections) {
  
  clusterIndicesUnique <- unique(clusterIndices[clusterIndices > 0])
  numClusters <- length(clusterIndicesUnique)
  
  clusterIndicesFilt <- rep(0, length(clusterIndices))
  
  clusterCount <- 0
  if(numClusters > 0) {
    for (i in 1 : numClusters) {
    
      numDetectionsCluster <- sum(clusterIndices == clusterIndicesUnique[i])
   
      if (numDetectionsCluster >= minDetections) {
        clusterCount <- clusterCount + 1
        clusterIndicesFilt[clusterIndices == clusterIndicesUnique[i]] <- clusterCount
      }
      
    }
  }
  return(clusterIndicesFilt)
}

