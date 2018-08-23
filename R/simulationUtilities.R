


#' Simulates dSTORM imaging process to generate synthetic single molecule localization microscopy datasets
#' 
#' @param moleculeCoords A matrix containing x and y (and z in 3D) coordinates of the underlying molecule distribution.
#' @param averageLabelsPerMol Average number of fluophores attatched to each molecule.
#' @param blinkGeomProb Probablilty of transistion to a dark state. Used to generate a geometric distrubtion.
#' @param precisionMeanLog Mean for lognormal distrubtion used to determine molecule localisation precision.
#' @param precisionSdLog Standard deviation for lognormal distrubtion used to determine molecule localisation precision.
#' @param detectionRate Percentage of blinking events which are detected.
#' @param falseDetectionRate Percentage of noise detections to add.
#' @param fieldLimit A matrix containing the field size limits.
#' @return Data frame containing the detection list (coordinates for all detections), and also the molecule list (linking detections to the orginating molecules)



simulateSTORM <- function(moleculeCoords, averageLabelsPerMol, blinkGeomProb, precisionMeanLog, precisionSdLog, detectionRate, falseDetectionRate, fieldLimit) {
  
  # number of molecules coordinates provided
  numMolecules <- dim(moleculeCoords)[1]
  
  # dimensionallity of the data
  numDimensions <- dim(moleculeCoords)[2]
  
  if (numDimensions > 3 || numDimensions < 2) {
    stop('Molecule coordinates should be 2D or 3D')
  } 
  
  # number of flurophores to be disitrubted between molecules
  numLabels <- averageLabelsPerMol * numMolecules
  
  # labels (flurophores) assigned to moelcules randomly
  labelIndex <- floor(runif(numLabels, min=1, max=numMolecules + 1))
  
  # retrieve label coordinates
  labelCoords <- matrix(0, nrow=numLabels, ncol=numDimensions)
  for (i in 1 : numLabels) {
    for (d in 1 : numDimensions) {
      labelCoords[i, d] <- moleculeCoords[labelIndex[i], d]
    }
  }

  # use a geometric distribution to determine the number of blinking events per label
  numBlinksPerLab <- rgeom(numLabels, blinkGeomProb)
  
  # total number of blinking events
  numBlinks = sum(numBlinksPerLab)
  # creat detection lists
  
  # to hold the molecule reference for each detection
  moleculeIndex <- c()
  # to hold the coordinates for detections
  detectionCoords <- c()
  
  for (i in 1 : numLabels) {
    
    # add the molecule indexes for this label
    moleculeIndex <- c(moleculeIndex, rep(labelIndex[i], numBlinksPerLab[i]))
    
    # to hold detection coordiantes for this label
    detectionCoordsLab = matrix(0, nrow=numBlinksPerLab[i], ncol=numDimensions)
    for (d in 1 : numDimensions) {
      # determine standard deviation for localization precision using lognormal distrubtion
      locPrec <- rlnorm(1, meanlog = precisionMeanLog[d], sdlog = precisionSdLog[d])
      # determine detection coordiantes using label coordinates and normal distribution 
      detectionCoordsLab[ , d] <- rep(labelCoords[i, d], numBlinksPerLab[i]) + rnorm(numBlinksPerLab[i], mean = 0, sd = locPrec)
    }
    
    # add to detection list
    detectionCoords <- rbind(detectionCoords, detectionCoordsLab)  
  }
  
  
  # remove missed detections
  detectionsIndicesKeep <- sample(1:numBlinks, round(detectionRate * numBlinks), replace=FALSE) 
  detectionCoords <- detectionCoords[detectionsIndicesKeep, ]
  moleculeIndex <- moleculeIndex[detectionsIndicesKeep]
  
  # remove detections which now lie outside the field of view
  detectionsInside = rep(TRUE, length(moleculeIndex))
  for (d in 1 : numDimensions) {
    detectionsInside[detectionCoords[ , d] < fieldLimits[1, d]] <- FALSE
    detectionsInside[detectionCoords[ , d] > fieldLimits[2, d]] <- FALSE
  }
  detectionCoords <- detectionCoords[detectionsInside, ]
  moleculeIndex <- moleculeIndex[detectionsInside]
  
  # add false detections (noise)
  numFlaseDetections <- round(falseDetectionRate * length(moleculeIndex))
  
  if (numFlaseDetections > 0) {
    # to hold coordinates of false detections
    falseDetections <- matrix(0, nrow=numFlaseDetections, ncol=numDimensions)
    for (d in 1 : numDimensions) {
      # place randomly on field of view
      falseDetections[, d] <-  runif(numFlaseDetections, fieldLimits[1, d], fieldLimits[2, d])
    }
    # assigned molecule index value of 0 and add to detection list
    moleculeIndex <- c(moleculeIndex, rep(0, numFlaseDetections))
    detectionCoords <- rbind(detectionCoords, falseDetections)
    
  }
  
  # build and return output data frame
  if (numDimensions == 2) {
    detectionList <- data.frame(x = detectionCoords[ , 1], y = detectionCoords[ , 2], moleculeIndex = moleculeIndex)
  } else {
    detectionList <- data.frame(x = detectionCoords[ , 1], y = detectionCoords[ , 2], z =detectionCoords[ , 3], moleculeIndex = moleculeIndex)
  }
  
  return (detectionList)
}


