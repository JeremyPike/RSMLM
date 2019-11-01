library("TDA")
library("plyr")

#' Utility function to return the number of features above a specified threshold from a persistence diagram.
#'
#' @param diag 3 Column matrix containing the perisitence diagram.  Rows correspond to individual features. Columns correspond to feature dimension, birth scale and death scale.
#' @param dimension The dimension of features to search for (0: connected components, 1: holes, 2: voids).
#' @param threshold Persitence threshold for features.
#' @return The number of features of the specified dimension above the specified threshold.
#' 

findNumFeatures <- function(diag, dimension, threshold) {
 
  if(dim(diag)[2] != 3) {
    stop('Diagram should have three columns')
  }

  numFilteredFeatures <- 0

  
  # remove features from the diagram which are not of the specified dimension
  diagDim <- diag[diag[ ,1]==dimension,]
  
  # if there is only one feature check if it is above the threshold
  if (is.vector(diagDim)) {
    numFeatures <- 1
    # persistence is birth scale - diagDim scale
    featurePersistence <- diagDim[3] - diag[2]
    if (featurePersistence > threshold) {
      numFilteredFeatures <- 1
    }
  } else {
    numFeatures <- dim(diagDim)[1]
    
    if (numFeatures != 0) {
      for (i in 1 : numFeatures) {
        # persistence is birth scale - death scale
        featurePersistence <- diagDim[i, 3] - diagDim[i, 2]
        # if persitence above threshold add one to feature count
        if (featurePersistence > threshold) {
          numFilteredFeatures <- numFilteredFeatures + 1
        }
      }
    }
  }
  
  return(numFilteredFeatures)
}

  


#' Sample detection coordinates (with) replacement a set number of times. For each set of sample coordinates compute the persitence diagram.

#' @param coords A matrix containing coordinates of the detections.
#' @param maxDimension The maximum dimenstion of topological features to calculate (0: CCs, 1: CCs/holes, 2: CCs/holes/voids).
#' @param maxScale Maximum scale over which to calculate the filtration
#' @param numBootstraps number of repeats for sub-sampling process. Default 20.
#' @param samplingRate A number between 0 and 1 specify the rate at which to sample the coordinates. Default 1.
#' @param weightingParameter (Optional) Parameter used to weight the sampling. For example the localization precision.
#' @param weightingFactor Posive number specifying the weighting factor for exponential decay function. Default log(0.1).
#' @param numberInGroup (Optional) Number of detections grouped into each specified coordinate.
#' @return List containing the persistent diagrams from the sampled data.
#' 
#' 
sampledDiagrams <- function(coords, maxDimension, maxScale, numBootstraps = 50, samplingRate = 1, weightingParameter, weightingFactor = 2.302585, numberInGroup) {
  
  
  
  numDimensions <- dim(coords)[2]
  if (numDimensions > 3 || numDimensions < 2) 
    stop('Coordinates should be 2D or 3D')
  
  if (maxDimension + 1 > numDimensions  || maxDimension %% 1 != 0 || maxDimension <= 0)
    stop('maxDimension should be a positive integer not exceeding dimensionality of data')
  
  if(numBootstraps %% 1 != 0 || numBootstraps <= 0)
    stop('numBootstraps should be an positive integer')
  
  if(samplingRate > 1 || samplingRate < 0)
    stop('sampling rate should be a value between 0 and 1')
  
  if(weightingFactor <= 0) 
    stop('weightingFactor whould be positive')
  
  if(!missing(numberInGroup)) {
    if(length(numberInGroup) != dim(coords)[1] || !all(numberInGroup == floor(numberInGroup)) || !all(numberInGroup > 0))
      stop('numberInGroup should be a vector with same number of elements as rows in coords. Each element should be a postive integer')
    # replicate coordinates and corresponding weighting parameter when detections have been grouped
    groupIndices <- unlist(mapply(FUN = rep, seq_along(numberInGroup), numberInGroup))
    coords <- coords[groupIndices, ]
    if (!missing(weightingParameter))
      weightingParameter <- weightingParameter[groupIndices]
  } 
  
  numDetections <- dim(coords)[1]
  
  # calculate the number of coordinates to sample
  m <- ceiling(numDetections / samplingRate)
  
  if(missing(weightingParameter)) {
    # if no weighting paramter is specified
    w <- rep(1, numDetections)
  } else {
    # normalise weighting paramter between 0 and 1
    normParameter <- (weightingParameter - min(weightingParameter)) / (max(weightingParameter) - min(weightingParameter))
    # calculate weighting using negative exponetial of normalised parameters
    w <- exp(-weightingFactor * normParameter)
  }

  diags <- list()
  for (j in seq_len(numBootstraps)) { 
    
    #print(paste("bootstrap", j, "of" , numBootstraps))
    # sample coordinates with replacement and specified weighting
    sampledCoords <- unique(coords[sample(seq_len(numDetections), m, replace = TRUE, prob = w), ])
    # calculate rips diagram for sampled data
    diags [[j]] <- ripsDiag(sampledCoords, maxdimension = maxDimension, maxscale = maxScale)$diagram 
    
  }
  

  return(diags)
}

#' Givin a list of diagrams calculated the number of features above a specified threshold and calculates the modal value and frequency
#'
#' @param diag 3 Column matrix containing the perisitence diagram.  Rows correspond to individual features. Columns correspond to feature dimension, birth scale and death scale.
#' @param dimension The dimension of features to search for (0: connected components, 1: holes, 2: voids).
#' @param threshold Persitence threshold for features.
#' @return List containing the modal number of features (numFeatures) and the fraction of diagrams which find the modal number (agreement)
#' 

findNumFeaturesConsensus <- function(diags, dimensions, thresholds) {
  
  numDimensions <- length(dimensions)
  if(numDimensions != length(thresholds)) {
    stop('The same number of dimensions and thresholds should be provided')
  }
  
  numDiagrams <- length(diags)
  
  numFeatures <- ldply(diags, .fun = function(x) {
    sapply(seq(1, numDimensions), FUN = function(y) {
      findNumFeatures(x, dimension = dimensions[y], threshold = thresholds[y])
    })
  })
  
  configCounts <- count(numFeatures)
  
  configConsenus <- configCounts[which(configCounts$freq == max(configCounts$freq))[1], ]
  
  agreement <- configConsenus$freq / numDiagrams
  
  return(list(numFeatures = as.numeric(configConsenus[1:numDimensions]), agreement = agreement))
}


