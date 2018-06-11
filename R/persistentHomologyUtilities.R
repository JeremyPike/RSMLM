
PH <- function(coords, maxDimension, maxScale, maxBottleneckDimension, numBootstraps) {
  
  #calculate persistence diagram
  diag <- ripsDiag(X = coords, maxDimension, maxScale, library = "GUDHI", printProgress = TRUE, location=TRUE)
  numDetections <- dim(coords)[1]
  m <- ceiling(numDetections / log(numDetections))
  bottleneckDist <- matrix(nrow = numBootstraps, ncol = maxDimension + 1)
  
  for (j in seq_len(numBootstraps)) { 
    
    print(paste("bootstrap", j, "of" , numBootstraps))
    
    # sample data
    subX <- coords[sample(seq_len(numDetections), m), ] 
    # calculate rips diagram for sampled data
    diagSamp <- ripsDiag(subX, maxdimension = maxDimension, maxscale = maxScale) 
    
    # calcualte bottneck distance for first two dimensions and also combined
    
    for (k in 0 : maxBottleneckDimension) {

      bottleneckDist[j, k + 1] <- bottleneck(diag$diagram, diagSamp$diagram, dimension = k)
    }
  }
  output <- list(diag = diag$diagram, bottleneckDist = bottleneckDist)
  return(output)
}



findNumFeatures <- function(diag, dimension, threshold) {
 
  numFilteredFeatures <- 0
  diag <- diag[diag[ ,1]==dimension,]
  if (is.vector(diag)) {
    numFeatures <- 1
  } else {
    numFeatures <- dim(diag)[1]
  }
  if (numFeatures != 0) {
    for (h in 1 : numFeatures) {
      if (is.vector(diag)) {
        featurePersistence <- diag[3] - diag[2]
      } else {
        featurePersistence <- diag[h, 3] - diag[h, 2]
      }
      if (featurePersistence > threshold) {
        numFilteredFeatures <- numFilteredFeatures + 1
      }
    }
  }
  return(numFilteredFeatures)
}