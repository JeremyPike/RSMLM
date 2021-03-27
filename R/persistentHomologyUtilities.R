

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
  
  # count to hold the number of filtered features
  numFilteredFeatures <- 0
  
  # remove features from the diagram which are not of the specified dimension
  diag <- diag[diag[ ,1]==dimension,]
  
  # if there is only one feature check if it is above the threshold
  if (is.vector(diag)) {
    numFeatures <- 1
    # persistence is birth scale - death scale
    featurePersistence <- diag[3] - diag[2]
    if (featurePersistence > threshold) {
      numFilteredFeatures <- 1
    }
  } else {
    numFeatures <- dim(diag)[1]
  }
  
  if (numFeatures != 0) {
    if (numFeatures > 1){
      for (i in 1 : numFeatures) {
        # persistence is birth scale - death scale
        featurePersistence <- diag[i,3] - diag[i,2]
        # if persitence above threshold add one to feature count
        if (featurePersistence > persistenceThreshold) {
          numFilteredFeatures <- numFilteredFeatures + 1
        }
      }
    }
    else{
      # persistence is birth scale - death scale
      featurePersistence <- diag[3]-diag[2]
      # if persitence above threshold add one to feature count
      if (featurePersistence > persistenceThreshold) {
        numFilteredFeatures <- numFilteredFeatures + 1
      }
    }
    
  }
  return(numFilteredFeatures)
}

