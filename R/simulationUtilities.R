


#' Simulates dSTORM imaging process to generate synthetic single molecule localization microscopy datasets
#' 
#' @param moleculeCoords A matrix containing x and y (and z in 3D) coordinates of the underlying molecule distribution.
#' @param averageLabelsPerMol Average number of fluophores attatched to each molecule.
#' @param blinkGeomProb Probablilty of transistion to a dark state. Used to generate a geometric distrubtion.
#' @param precisionMeanLog Mean for lognormal distrubtion used to determine molecule localisation precision.
#' @param precisionSdLog Standard deviation for lognormal distrubtion used to determine molecule localisation precision.
#' @param detectionRate Percentage of blinking events which are detected.
#' @param falseDetectionRate Percentage of noise detections to add.
#' @param fieldLimits A matrix containing the field size limits.
#' @param precisionMeanLogNoise (Optional) Mean for lognormal distrubtion used to determine localisation precision of noisy detections.
#' @param precisionSdLogNoise (Optional) Standard for lognormal distrubtion used to determine localisation precision of noisy detections.
#' @return Data frame containing the detection list (coordinates for all detections), and also the molecule list (linking detections to the orginating molecules)


simulateSTORM <- function(moleculeCoords, averageLabelsPerMol, blinkGeomProb, precisionMeanLog, precisionSdLog, detectionRate, falseDetectionRate, fieldLimits, precisionMeanLogNoise, precisionSdLogNoise) {
  
  
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
  labelIndex <- sample(1:numMolecules, numLabels,  TRUE)
  
  # retrieve label coordinates
  labelCoords <- matrix(0, nrow=numLabels, ncol=numDimensions)
  for (i in 1 : numLabels) {
    for (d in 1 : numDimensions) {
      labelCoords[i, d] <- moleculeCoords[labelIndex[i], d]
    }
  }
  
  detectionList <- simulateSMLM(labelCoords, blinkGeomProb, precisionMeanLog, precisionSdLog, detectionRate, falseDetectionRate, fieldLimits, precisionMeanLogNoise, precisionSdLogNoise, labelIndex)
 
  
  return (detectionList)
}


#' Simulates SMLM imaging process to generate synthetic single molecule localization microscopy datasets.
#' 
#' @param labelCoords A matrix containing x and y (and z in 3D) coordinates of the underlying label distribution.
#' @param blinkGeomProb Probablilty of transistion to a dark state. Used to generate a geometric distrubtion.
#' @param precisionMeanLog Mean for lognormal distrubtion used to determine molecule localisation precision.
#' @param precisionSdLog Standard deviation for lognormal distrubtion used to determine molecule localisation precision.
#' @param detectionRate Percentage of blinking events which are detected.
#' @param falseDetectionRate Percentage of noise detections to add.
#' @param fieldLimits A matrix containing the field size limits.
#' @param precisionMeanLogNoise (Optional) Mean for lognormal distrubtion used to determine localisation precision of noisy detections.
#' @param precisionSdLogNoise (Optional) Standard for lognormal distrubtion used to determine localisation precision of noisy detections.
#' @param labelIndex (Optional) Vector linking labels to molecules
#' @return Data frame containing the detection list (coordinates for all detections), and also the molecule list (linking detections to the orginating molecules)



simulateSMLM <- function(labelCoords, blinkGeomProb, precisionMeanLog, precisionSdLog, detectionRate, falseDetectionRate, fieldLimits, precisionMeanLogNoise, precisionSdLogNoise, labelIndex) {
  
  # number of label coordinates provided
  numLabels <- dim(labelCoords)[1]
  
  # dimensionallity of the data
  numDimensions <- dim(labelCoords)[2]
  
  if(missing(precisionMeanLogNoise))
    precisionMeanLogNoise <- precisionMeanLog
  if(missing(precisionSdLogNoise))
    precisionSdLogNoise <- precisionSdLog
  if(missing(labelIndex))
    labelIndex <- 1:numLabels

  if (numDimensions > 3 || numDimensions < 2) {
    stop('Molecule coordinates should be 2D or 3D')
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
  
  # to hold errors (absolute value of localization precison)
  errors <- c()
  
  for (i in 1 : numLabels) {
    
    # add the molecule indexes for this label
    moleculeIndex <- c(moleculeIndex, rep(labelIndex[i], numBlinksPerLab[i]))
    
    # to hold detection coordiantes for this label
    detectionCoordsLab = matrix(0, nrow=numBlinksPerLab[i], ncol=numDimensions)
    locPrec = matrix(0, nrow=numBlinksPerLab[i], ncol=numDimensions)
    for (d in 1 : numDimensions) {
      # determine standard deviation for localization precision using lognormal distrubtion
      locPrecSd <- rlnorm(1, meanlog = precisionMeanLog[d], sdlog = precisionSdLog[d])
      locPrec[ , d] <- rnorm(numBlinksPerLab[i], mean = 0, sd = locPrecSd)
      # determine detection coordiantes using label coordinates and normal distribution 
      detectionCoordsLab[ , d] <- rep(labelCoords[i, d], numBlinksPerLab[i]) + locPrec[ , d]
    }
    
    # add to detection list
    detectionCoords <- rbind(detectionCoords, detectionCoordsLab)  
    errors <- rbind(errors, abs(locPrec))
  }
  
  
  # remove missed detections
  detectionsIndicesKeep <- sample(1:numBlinks, round(detectionRate * numBlinks), replace=FALSE) 
  detectionCoords <- detectionCoords[detectionsIndicesKeep, ]
  errors <- errors[detectionsIndicesKeep, ]
  moleculeIndex <- moleculeIndex[detectionsIndicesKeep]
  
  # remove detections which now lie outside the field of view
  detectionsInside = rep(TRUE, length(moleculeIndex))
  for (d in 1 : numDimensions) {
    detectionsInside[detectionCoords[ , d] < fieldLimits[1, d]] <- FALSE
    detectionsInside[detectionCoords[ , d] > fieldLimits[2, d]] <- FALSE
  } 

  detectionCoords <- detectionCoords[detectionsInside, ]
  errors <- errors[detectionsInside, ]
  moleculeIndex <- moleculeIndex[detectionsInside]

  # add false detections (noise)
  numFlaseDetections <- round(falseDetectionRate * length(moleculeIndex))
  
  if (numFlaseDetections > 0) {
    # to hold coordinates of false detections
    falseDetections <- matrix(0, nrow=numFlaseDetections, ncol=numDimensions)
    locPrec = matrix(0, nrow=numFlaseDetections, ncol=numDimensions)
    for (d in 1 : numDimensions) {
      # place randomly on field of view
      falseDetections[, d] <-  runif(numFlaseDetections, fieldLimits[1, d], fieldLimits[2, d])
      locPrecSd <- rlnorm(numFlaseDetections, meanlog = precisionMeanLogNoise[d], sdlog = precisionSdLogNoise[d])
      locPrec[ ,d] <- rnorm(numFlaseDetections, mean = 0, sd = locPrecSd)
    }
    # assigned molecule index value of 0 and add to detection list
    moleculeIndex <- c(moleculeIndex, rep(0, numFlaseDetections))
    detectionCoords <- rbind(detectionCoords, falseDetections)
    errors <- rbind(errors, abs(locPrec))
  }
  

  errorsRMS <- sqrt(rowSums(errors^2))

  # build and return output data frame
  if (numDimensions == 2) {
    detectionList <- data.frame(x = detectionCoords[ , 1], y = detectionCoords[ , 2], moleculeIndex = moleculeIndex, error = errorsRMS)
  } else {
    detectionList <- data.frame(x = detectionCoords[ , 1], y = detectionCoords[ , 2], z =detectionCoords[ , 3], moleculeIndex = moleculeIndex, error = errorsRMS)
  }
  
  return (detectionList)
  
  
  
}