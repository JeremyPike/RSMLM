
library("dbscan")
library("deldir")
library("igraph")
library("geometry")
library("ggplot2")
library("randomcoloR")



#' Topological Model Analysis Tool (ToMATo) clustering
#'
#' As orginally described in Chazal, Frederic, et al. "Persistence-based clustering in riemannian manifolds." Journal of the ACM (JACM) 60.6 (2013): 41.
#' 
#' @param coords A matrix containing coordinates of the detections.
#' @param r The search radius for density estimation and to contruct linking graph
#' @param threshold Persistence threshold for merging of density modes
#' @return Vector containing the cluster indices for each detection, a cluster index of zero refers to noise
#' 
clusterTomato <- function(coords, r, threshold) {
  
  coords <- as.matrix(coords)
  
  numDimensions <- dim(coords)[2]
  if (numDimensions > 3 || numDimensions < 2) {
    stop('Coordinates should be 2D or 3D')
  } 
  
  # compute density using number of other detections within fixed search radius (library dbscan)
  density <- pointdensity(coords, r, type = "frequency")
  # tomato algorithm performed using c++ function
  clusterIndices <- tomatoDens(coords, density, r, threshold)$clusters
  return(clusterIndices)
}



#' Density-based spatial clustering of applications with noise (DBSCAN)
#' 
#' As orginally described in Ester, Martin, et al. "Density-based spatial clustering of applications with noise." Int. Conf. Knowledge Discovery and Data Mining. Vol. 240. 1996.
#' This function simply calls the dbscan function from the library dbscan and is included for convenience. Border detections are included in clusters. 
#' 
#' @param coords A matrix containing coordinates of the detections.
#' @param r The search radius for density estimation 
#' @param minPoints Minimum number of detections that should be within the search radius for a detection to be defined as a "core" point.
#' @return Vector containing the cluster indices for each detection, a cluster index of zero refers to noise
#' 
clusterDBSCAN <- function(coords, r, minPoints) {
  
  numDimensions <- dim(coords)[2]
  if (numDimensions > 3 || numDimensions < 2) {
    stop('Coordinates should be 2D or 3D')
  } 
  # perform DBSCAN using "dbscan" pacakage
  fr <- frNN(coords, r)
  clusterIndices <- dbscan(fr, minPts = minPoints, borderPoints=TRUE)$cluster
  
  return(clusterIndices)
}


#' 2D Voronoi tesselation based clustering
#' 
#' A Voronoi tesselation is contructed from the detection coordinates. Detection density is defined as the inverse of tile area, or mean first rank tile area.
#' Detections with density below a specified threshold are removed from the dual Delaunay triangulation and the connected components of this graph are used to define clusters. 
#' Implementation is for 2D data only. \cr \cr
#' Levet, Florian, et al. "SR-Tesseler: a method to segment and quantify localization-based super-resolution microscopy data." Nature methods 12.11 (2015): 1065. \cr
#' Andronov, Leonid, et al. "ClusterViSu, a method for clustering of protein complexes by Voronoi tessellation in super-resolution microscopy." Scientific reports 6 (2016): 24084.
#' 
#' @param coords A matrix containing coordinates of the detections.
#' @param threshold Threshold for density.
#' @param densityChoice Choice of density calculation method. 0: inverse tile area, 1: inverse mean first rank tile area
#' @return Vector containing the cluster indices for each detection, a cluster index of zero refers to noise
#'
clusterVoronoi <- function(coords, threshold, densityChoice) {
  
  numDimensions <- dim(coords)[2]
  if (numDimensions != 2) {
    stop('Coordinates should be 2D')
  } 
  if (sum(duplicatedxy(coords[, 1], coords[, 2])) > 0) {
    stop("Coords list should not contain duplicated points")
  } 
  
  numDetections <- dim(coords)[1]
  
  # compute veronoi tesselation using deldir library
  vtess <- deldir(coords[, 1], coords[, 2])
  
  # retrieve area of Voronoi tiles
  tileAreas <- vtess$summary$dir.area
 
  
  # create igraph object corresponding to the Delaunay triangulation 
  g <- make_empty_graph(n = numDetections, directed = FALSE)
  g <- add_edges(g, c(rbind(vtess$delsgs$ind1, vtess$delsgs$ind2)))
  
  # to hold density estimate
  density <- rep(0, numDetections)
  
  # calculate neighbours for each detection in the graph
  neigh <- adjacent_vertices(g, seq(1, numDetections, by = 1))
  
  # use first rank tile area for density calculation
  if (densityChoice == 1) {
    for (i in 1 : numDetections) {
      density[i] <- (1 + length(neigh[[i]])) / (tileAreas[i] + sum(tileAreas[neigh[[i]]]))
    }
  # use tile area for density calculation
  } else if (densityChoice == 0) {
    for (i in 1 : numDetections) {
      density[i] <- (1 + length(neigh[[i]])) / (tileAreas[i] + sum(tileAreas[neigh[[i]]]))
      
    }
  }
  # delete tiles with area larger than specified maxiumum
  tileIndices <-  density > threshold
  g <- delete_vertices(g, which(!tileIndices))
  
  # compute the connected components of the graph
  ccs <- components(g)
  # extract cluster indices as vector with zero value representing no cluster assignment
  clusterIndices <- rep(0, numDetections)
  clusterIndices[which(tileIndices)] <- ccs$membership
  
  return(clusterIndices)
  
}


#' Ripley's K function based clustering
#' 
#' First the number of detections within a fixed search radius, r, is counted. This is used to claculate the Ripley K function which is normalised to get the Ripley L function.
#' Detections with L-r less than a specified threshold are removed. Clusters are defined by the connected components of the graph formed by linking all filtered detections within r.
#' This is approximatelly as described in the following paper but without Baysian parameter selection: \cr \cr
#' Rubin-Delanchy, Patrick, et al. "Bayesian cluster identification in single-molecule localization microscopy data." Nature methods 12.11 (2015): 1072.
#' 
#' @param coords A matrix containing coordinates of the detections.
#' @param r The search radius for calcualte of the Ripley K and L functions.
#' @param threshold Threshold for L-r.
#' @param ROIArea Area of the field of view. Needed to calcualte the Ripley K function.
#' @return Vector containing the cluster indices for each detection, a cluster index of zero refers to noise
#'

clusterRipley <- function(coords, r, threshold, ROIArea) {
  
  numDimensions <- dim(coords)[2]
  if (numDimensions > 3 || numDimensions < 2) {
    stop('Coordinates should be 2D or 3D')
  } 
  
  numDetections <- dim(coords)[1]
  
  # use dbscan library to count number of neighbours within fixed search radius
  fr <- frNN(coords, eps = r)
  numNeighbours <- rep(0, numDetections)
  for (k in 1 : numDetections) {
    numNeighbours[k] <- length(fr$id[[k]])
  }
  
  # calculate the Ripley K function for each detection
  K <- numNeighbours * ROIArea / (numDetections - 1)
  
  # calculate the Ripley L function and subract r to get L - r (LminR)
  if (numDimensions == 3) {
    LminR <- (3 * K / pi / 4)^(1/3) - r
  } else {
    LminR <- sqrt(K / pi) - r
  }
  
  # to hold custer indices
  clusterIndices <- rep(0, numDetections)
  
  # find which detections are above the specified threshold
  filteredDetections <- LminR > threshold
  numFilteredDetections <- sum(filteredDetections)
  
  # if just one detction is filtered then assign this a cluster index of 1
  if(numFilteredDetections == 1) {

    clusterIndices[filteredDetections] <- 1
  } else if (numFilteredDetections > 1)  {
    # contruct graph using r to link detections
    g <- graph_from_adj_list((adjacencylist)(fr))
    # delete vertices with L-r below the specified threshold
    g <- delete_vertices(g, which(!filteredDetections)) 
    # find the connected compontents of the graph
    ccs <- components(g)
    
    # use connected components to define cluster indices
    clusterIndices[filteredDetections] <- ccs$membership
  }
  return(clusterIndices)
}

#' Simple clustering scheme linking all detections within a search radius
#' 
#' Forms a graph by linking all detections within a specified search radius. Clusters defined as connected components of this graph.
#' 
#' @param coords A matrix containing coordinates of the detections.
#' @param r The search radius used to contruct the graph.
#' @return Vector containing the cluster indices for each detection, a cluster index of zero refers to noise
#'
clusterCCs <- function(coords, r) {
  
  # fixed radius neighbour search
  fr <- frNN(coords, eps = r)
  # create graph
  g <- graph_from_adj_list(adjacencylist(fr))
  # find the connected components of the graph
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


#' Produces a ToMATo persistence diagram
#'
#' As orginally described in Chazal, Frederic, et al. "Persistence-based clustering in riemannian manifolds." Journal of the ACM (JACM) 60.6 (2013): 41.
#' Records the bith and death densities for all density modes
#' Calculated by running the ToMATo clustering with very high threshold
#' 
#' @param coords A matrix containing coordinates of the detections.
#' @param r The search radius for density estiamtion and to contruct linking graph
#' @return Matix containing the bith and death densities for all modes in the density estimate
#' 
tomatoDiagram <- function(coords, r) {
  
  coords <- as.matrix(coords)
  
  numDimensions <- dim(coords)[2]
  if (numDimensions > 3 || numDimensions < 2) {
    stop('Coordinates should be 2D or 3D')
  } 
  
  # compute density using number of other detections within fixed search radius (library dbscan)
  density <- pointdensity(coords, r, type = "frequency")
  # computes the digram using C++ scripts
  diagram <- tomatoDens(coords, density, r, 1e20)$diagram
  return(diagram)
}


#' Simple plotting function for ToMATo diagrams
#'
#' 
#' @param diagram Matix containing the bith and death densities for all modes in the density estimate
#' @param threshold Persistence threshold for merging of density modes
#' 
plotTomatoDiagram <- function(diagram, threshold) {

  # maxium birth density
  maxBirth <- max(diagram[ , 1])
  
  # plot all modes and limit axis using the maximum birth density
  plot(diagram[ , 1], diagram[ , 2], xlab="Birth density", ylab="Death density", asp=1, xlim=c(0, maxBirth), ylim=c(-1, maxBirth), pch = 16, cex = .8)
  # plot line through the diagonal
  lines(c(0, maxBirth), c(0, maxBirth))
  # plot location of specified persistence threshold
  lines(c(threshold, maxBirth), c(0, maxBirth - threshold), lty=2)
  lines(c(threshold, threshold), c(-1, 0), lty=2)
  
  
}


#' A simple 2D plotting function using ggplot2 to display clustering results.
#'
#' Distinct cluster colours are sampled from a selection produced by randomColoR
#' 
#' @param coords A matrix containing x and y coordinates of the detections.
#' @param clusterIndices The cluster index for each detection. Noise detections should have cluster index zero.


plotClusterScatter <- function(coords, clusterIndices) {
  
  # dimensionallity of the data
  numDimensions <- dim(coords)[2]

  if (numDimensions != 2) {
    stop('Data should be 2D')
  } 
  
  # check coords and clsterIndices have the same number of detections
  numDetections <- dim(coords)[1]
  if (numDetections != length(clusterIndices)) {
    stop('coords and clusterIndices should have the same length')
  }   
    
  # find unique clusters
  uniqueClusters <- unique(clusterIndices)

  numClusters <- sum(uniqueClusters > 0)
  
  # contruct data frame from input data  
  detectionList <- data.frame(x = coords[ , 1], y = coords[ , 2], clusterIndex = clusterIndices)
  
  # if there are clusters
  if (numClusters > 0) {
    
    # generate cluster colours
    colorMap <- as.character(distinctColorPalette(numClusters), luminosity = "dark")
  
    #if there are noise detections
    if (0 %in% uniqueClusters) {
      # plot these in black and ranomise order of other clusters
      clusterColors <- c("black", sample(colorMap, numClusters))
    } else {
      clusterColors <- sample(colorMap, numClusters)
    }
    
    # plot using ggplot2
    ggplot(detectionList, aes(x = x, y = y, color = factor(clusterIndex))) + 
      geom_point(size = 1.2, shape = 16) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = clusterColors) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            legend.position="none")  
  } else {
    ggplot(detectionList, aes(x = x, y = y)) + 
      geom_point(size = 1.2, shape = 16, color = "black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            legend.position="none") 
    
  }
  
}



#' Calculate per cluster statsitics
#'
#' Calculates the number of detections, area, volume density and centre of mass for each cluster. Returns statistics for each cluster and also mean values.
#' 
#' @param coords A matrix containing coordinates of the detections.
#' @param clusterIndices The cluster index for each detection. Noise detections should have cluster index zero.
#' @return A data frame containing statistics for each cluster
#' 


clusterStats <- function(coords, clusterIndices) {
  
  numDimensions <- dim(coords)[2]
  if (numDimensions > 3 || numDimensions < 2) {
    stop('Coordinates should be 2D or 3D')
  } 
  
  # check coords and clsterIndices have the same number of detections
  numDetections <- dim(coords)[1]
  if (numDetections != length(clusterIndices)) {
    stop('coords and clusterIndices should have the same length')
  }   
  
  
  # find indices of all clusters
  clusterIndicesUnique <- unique(clusterIndices[clusterIndices > 0])
  numClusters <- length(clusterIndicesUnique)
  
  # data frame to hold per cluster statistics (fill with zeros)
  clustStats <-data.frame(matrix(0, nrow = numClusters, ncol = 7))
  colnames(clustStats) <- c("numDetectionsCluster", "areasCluster", "volumesCluster", "densitiesCluster", "meanX", "meanY", "meanZ")
  
  if (numClusters > 0) {

    for(i in 1 : numClusters) {
      
      # cluster coordinates
      coordsCluster <- coords[clusterIndices == i, ]
      
      # mean coordinates for each dimension
      clustStats$meanX[i]  <- mean(coordsCluster[, 1]) 
      clustStats$meanY[i]  <- mean(coordsCluster[, 2])
      if (numDimensions == 3) {
        clustStats$meanZ[i]  <- mean(coordsCluster[, 3])
      }
      
      # find number of detections in cluster
      clustStats$numDetectionsCluster[i] <- sum(clusterIndices == clusterIndicesUnique[i])
      
      
      if (clustStats$numDetectionsCluster[i] > 2) {
        #calculate convex hull
        ch <- convhulln(coordsCluster, options = "FA")
        # get cluster area (and also volume for 3D) and densities. Used geometry library to calculate convex hulls.
        if (numDimensions == 2) {
          clustStats$areasCluster[i] <- ch$vol
          clustStats$densitiesCluster[i] <- clustStats$numDetectionsCluster[i] / clustStats$areasCluster[i]
        } else {
          clustStats$areasCluster[i] <- ch$area
          clustStats$volumesCluster[i] <- ch$vol
          clustStats$densitiesCluster[i] <- clustStats$numDetectionsCluster[i] / clustStats$volumesCluster[i]
        }
        
      }

    }
  }

  return(clustStats)
}




#' Filter clusters based on number of detecions in a cluster
#'
#' Detections in clusters outside specified range are asigned an index of zero
#' 
#' @param clusterIndices The cluster index for each detection. Noise detections should have cluster index zero.
#' @param minDetections The minimum number of detections each cluster should contain.
#' @param maxDetections The maximum number of detections each cluster should contain.
#' @return Vector containing filtered cluster indices
#' 
filtClustDetections <- function(clusterIndices, minDetections, maxDetections) {
  
  # find indices of all clusters
  clusterIndicesUnique <- unique(clusterIndices[clusterIndices > 0])
  numClusters <- length(clusterIndicesUnique)
  
  # to hold filtered cluster indices
  clusterIndicesFilt <- rep(0, length(clusterIndices))
  
  # count of number of filtered clusters
  clusterCount <- 0
  
  if(numClusters > 0) {
    for (i in 1 : numClusters) {
      
      # find number of detections in cluster
      numDetectionsCluster <- sum(clusterIndices == clusterIndicesUnique[i])
      
      # if this in in the specified range
      if (numDetectionsCluster >= minDetections && numDetectionsCluster <= maxDetections) {
        clusterCount <- clusterCount + 1
        # use cluster count to define index for filtered cluster
        clusterIndicesFilt[clusterIndices == clusterIndicesUnique[i]] <- clusterCount
      }
      
    }
  }
  return(clusterIndicesFilt)
}

#' Filter clusters based on area/volume
#'
#'In 3D clusters are filtered by volume instead of area.
#'Area/volume defined by convex hull. 
#'Detections in clusters outside specified range are asigned an index of zero.
#' 
#' @param coords A matrix containing coordinates of the detections.
#' @param clusterIndices The cluster index for each detection. Noise detections should have cluster index zero.
#' @param minArea Minimum cluster area/volume.
#' @param maxArea Maximum cluster area/volume.
#' @return Vector containing filtered cluster indices
#' 
filtClustArea <- function(coords, clusterIndices, minArea, maxArea) {
  
  
  numDimensions <- dim(coords)[2]
  if (numDimensions > 3 || numDimensions < 2) {
    stop('Coordinates should be 2D or 3D')
  } 
  
  # check coords and clsterIndices have the same number of detections
  numDetections <- dim(coords)[1]
  if (numDetections != length(clusterIndices)) {
    stop('coords and clusterIndices should have the same length')
  }   
  
  # find indices of all clusters
  clusterIndicesUnique <- unique(clusterIndices[clusterIndices > 0])
  numClusters <- length(clusterIndicesUnique)
  
  # to hold filtered cluster indices
  clusterIndicesFilt <- rep(0, length(clusterIndices))
  
  # count of number of filtered clusters
  clusterCount <- 0
  
  if(numClusters > 0) {
    for (i in 1 : numClusters) {
      
      # find number of detections in cluster
      numDetectionsCluster <- sum(clusterIndices == clusterIndicesUnique[i])
      if (numDetectionsCluster > 2) {
        
        coordsCluster <- coords[clusterIndices == i, ]
        # calculate convex hull of cluster
        ch <- convhulln(coordsCluster, options = "FA")
        # for 2D data this is the area, for 3D the volume
        areaCluster <- ch$vol
        # if inside specified range
        if (areaCluster >= minArea && areaCluster <= maxArea) {
          clusterCount <- clusterCount + 1
          # use cluster count to define index for filtered cluster
          clusterIndicesFilt[clusterIndices == clusterIndicesUnique[i]] <- clusterCount
        }
      }
    }
  }
  return(clusterIndicesFilt)
}

