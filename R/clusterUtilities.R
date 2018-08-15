
library("dbscan")
library("deldir")
library("igraph")
library("geometry")

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
clusterVoronoi <- function(coords, threshold, algorithm, shortestDistance) {
  
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
    if (algorithm == 1) {
      density <- c()
      distances <- c()
      neigh <- adjacent_vertices(g, seq(1, numDetections, by = 1))
      for (i in 1 : numDetections) {
        density[i] <- (1 + length(neigh[[i]])) / (tileAreas[i] + sum(tileAreas[neigh[[i]]]))
        distances[i] <- min(sqrt((coords[i, 1] - coords[neigh[[i]], 1])^2 + (coords[i, 2] - coords[neigh[[i]], 2])^2))
      }
      tileIndices <-  density > threshold & distances < shortestDistance
    } else if (algorithm == 2) {
      density <- c()
      neigh <- adjacent_vertices(g, seq(1, numDetections, by = 1))
      for (i in 1 : numDetections) {
        density[i] <- (1 + length(neigh[[i]])) / (tileAreas[i] + sum(tileAreas[neigh[[i]]]))
        
      }
      tileIndices <- density > threshold
    }
    else {
      density <- 1 / tileAreas
      # find tiles with density larger than threshold
      tileIndices <- density > threshold
    }
 
    
    
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

plotTomatoDiagram <- function(diagram, threshold, main) {
  birthLim <- max(diagram[ , 1])
  plot(diagram[ , 1], diagram[ , 2], xlim=c(0, birthLim), ylim=c(-1, birthLim), type = "p", col = "black", pch = 19, cex = 0.5, xlab="birth", ylab="death", main=main)
  lines(c(0, birthLim),c(0, birthLim), type = "l", col="black")
  lines(c(threshold, birthLim),c(0, birthLim - threshold), type = "l", col="green")
  lines(c(threshold, threshold),c(-1, 0), type = "l", col="green")
  
  
}

filtClustDetections <- function(clusterIndices, minDetections, maxDetections) {
  
  clusterIndicesUnique <- unique(clusterIndices[clusterIndices > 0])
  numClusters <- length(clusterIndicesUnique)
  
  clusterIndicesFilt <- rep(0, length(clusterIndices))
  
  clusterCount <- 0
  if(numClusters > 0) {
    for (i in 1 : numClusters) {
    
      numDetectionsCluster <- sum(clusterIndices == clusterIndicesUnique[i])
   
      if (numDetectionsCluster >= minDetections && numDetectionsCluster <= maxDetections) {
        clusterCount <- clusterCount + 1
        clusterIndicesFilt[clusterIndices == clusterIndicesUnique[i]] <- clusterCount
      }
      
    }
  }
  return(clusterIndicesFilt)
}


filtClustAreas <- function(coords, clusterIndices, minArea, maxArea) {
  
  clusterIndicesUnique <- unique(clusterIndices[clusterIndices > 0])
  numClusters <- length(clusterIndicesUnique)
  
  clusterIndicesFilt <- rep(0, length(clusterIndices))
  clusterCount <- 0
  if(numClusters > 0) {
    for (i in 1 : numClusters) {
      
      numDetectionsCluster <- sum(clusterIndices == clusterIndicesUnique[i])
      if (numDetectionsCluster > 2) {
        coordsCluster <- coords[clusterIndices == i, ]
        ch <- convhulln(coordsCluster, options = "FA")
        areaCluster <- ch$vol
        if (areaCluster >= minArea && areaCluster <= maxArea) {
          clusterCount <- clusterCount + 1
          clusterIndicesFilt[clusterIndices == clusterIndicesUnique[i]] <- clusterCount
        }
      }
    }
  }
  return(clusterIndicesFilt)
}

plotClusterScatter <- function(coords, clusterIndices) {
  
  detectionCol <- rep(rgb(0,0, 0, maxColorValue = 255), dim(coords)[1])
  
  numClusters <- length(unique(clusterIndices)) - 1
  
  rbPal <- colorRampPalette(c('red','green'))
  
  if (numClusters > 0) {
    
    clusterCols <- sample(rbPal(numClusters))
    for (j in 1 : numClusters) {
      detectionCol[clusterIndices == j] <- clusterCols[j]
    }

  }
  plot(coords[ , 1], coords[ , 2], col = detectionCol, pch = 19)
}


clusterStats <- function(coords, clusterIndices) {

  numClusters <- length(unique(clusterIndices[clusterIndices>0]))
  numDetectionsCluster <- c()
  volumesCluster <- rep(0, numClusters)
  areasCluster <- rep(0, numClusters)
  densitiesCluster <- rep(0, numClusters)
  meanX  <- rep(0, numClusters)
  meanY  <- rep(0, numClusters)
  meanZ  <- rep(0, numClusters)

  if (numClusters > 0) {

    for(i in 1 : numClusters) {
      coordsCluster <- coords[clusterIndices == i, ]
      meanX[i]  <- mean(coordsCluster[, 1]) 
      meanY[i]  <- mean(coordsCluster[, 2])
  
      
      if (is.vector(coordsCluster)) {
        numDetectionsCluster[i] <- 1
      } else {
        numDetectionsCluster[i] <- dim(coordsCluster)[1]
      }
      if (numDetectionsCluster[i] > 2) {
        ch <- convhulln(coordsCluster, options = "FA")
        if (dim(coords)[2] == 2) {
          areasCluster[i] <- ch$vol
        } else {
          areasCluster[i] <- ch$area
          volumesCluster[i] <- ch$vol
          meanZ[i]  <- mean(coordsCluster[, 3])
        }
        
        densitiesCluster[i] <- numDetectionsCluster[i] / areasCluster[i] * 1000 * 1000
      }


    }
  }
  meanNumDetectionsCluster <- mean(numDetectionsCluster)
  sdNumDetectionsCluster <- sd(numDetectionsCluster)
  meanVolumesCluster <- mean(volumesCluster)
  sdVolumesCluster <- sd(volumesCluster)
  meanAreasCluster <- mean(areasCluster)
  sdAreasCluster <- sd(areasCluster)
  meanDensitiesCluster <- mean(densitiesCluster)
  sdDensitiesCluster <- sd(densitiesCluster)

  meanStats <- data.frame(numClusters, meanNumDetectionsCluster, sdNumDetectionsCluster,meanVolumesCluster, sdVolumesCluster, meanAreasCluster, sdAreasCluster, meanDensitiesCluster, sdDensitiesCluster)
  individualStats <- data.frame(numDetectionsCluster, areasCluster, volumesCluster, densitiesCluster, meanX, meanY, meanZ)
  stats <- list(meanStats = meanStats, individualStats = individualStats)

  return(stats)
}

