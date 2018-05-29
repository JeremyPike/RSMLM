simulateSTORM <- function(moleculeCoords, averageLabelsPerMol, blinkGeomProb, precisionMeanLog, precisionSdLog, detectionRate, falseDetectionRate, fieldLimitsX, fieldLimitsY) {
  
  numMolecules <- dim(moleculeCoords)[1]
  moleculePositionsX <- moleculeCoords[ , 1]
  moleculePositionsY <- moleculeCoords[ , 2]
  
  numLabelledMolecules <- averageLabelsPerMol * numMolecules
  
  labelIndex <- floor(runif(numLabelledMolecules, min=1, max=numMolecules + 1))
  labelPositionsX <- c()
  labelPositionsY <- c()

  
  for (i in 1 : numLabelledMolecules) {
    labelPositionsX[i] <- moleculePositionsX[labelIndex[i]]
    labelPositionsY[i] <- moleculePositionsY[labelIndex[i]]
  }
  
  numDetectionsPerMol <- rgeom(numLabelledMolecules, blinkGeomProb)
  # creat detection lists
  detectionPositionsX <- c()
  detectionPositionsY <- c()
  detectionIndex <- c()
  locPrec <- rlnorm(numLabelledMolecules, meanlog = precisionMeanLog, sdlog = precisionSdLog)
  for (i in 1 : numLabelledMolecules) {
    detectionPositionsXtemp <- rep(labelPositionsX[i], numDetectionsPerMol[i]) + rnorm(numDetectionsPerMol[i], mean = 0, sd = locPrec[i])
    detectionPositionsYtemp <- rep(labelPositionsY[i], numDetectionsPerMol[i]) + rnorm(numDetectionsPerMol[i], mean = 0, sd = locPrec[i])
    detectionPositionsX <- c(detectionPositionsX, detectionPositionsXtemp)  
    detectionPositionsY <- c(detectionPositionsY, detectionPositionsYtemp)
    detectionIndex <- c(detectionIndex, rep(labelIndex[i], numDetectionsPerMol[i]))
  }
  
 
  # remove missed detections
  detectionsIndices <- sample(1:length(detectionPositionsX), round(detectionRate * length(detectionPositionsX)), replace=FALSE) 
  detectionPositionsX <- detectionPositionsX[detectionsIndices]
  detectionPositionsY <- detectionPositionsY[detectionsIndices]
  detectionIndex <- detectionIndex[detectionsIndices]
  
  # remove detections which now lie outside the field of view
  detectionPositionsXFilt <- detectionPositionsX[detectionPositionsX >= fieldLimitsX[1] & detectionPositionsY >= fieldLimitsY[1]]
  detectionPositionsYFilt <- detectionPositionsY[detectionPositionsX >= fieldLimitsX[1] & detectionPositionsY >= fieldLimitsY[1]]
  detectionIndex <- detectionIndex[detectionPositionsX >= fieldLimitsX[1] & detectionPositionsY >= fieldLimitsY[1]]
  detectionPositionsX <- detectionPositionsXFilt
  detectionPositionsY <- detectionPositionsYFilt
  detectionPositionsXFilt <- detectionPositionsX[detectionPositionsX < fieldLimitsX[2] & detectionPositionsY < fieldLimitsY[2]]
  detectionPositionsYFilt <- detectionPositionsY[detectionPositionsX < fieldLimitsX[2] & detectionPositionsY < fieldLimitsY[2]]
  detectionIndex <- detectionIndex[detectionPositionsX < fieldLimitsX[2] & detectionPositionsY < fieldLimitsY[2]]
  detectionPositionsX <- detectionPositionsXFilt
  detectionPositionsY <- detectionPositionsYFilt
  
  # add false detections
  numFlaseDetections <- round(falseDetectionRate * length(detectionPositionsX))
  print(numFlaseDetections)
  if (numFlaseDetections > 0) {
    detectionPositionsX <- c(detectionPositionsX, runif(numFlaseDetections, fieldLimitsX[1], fieldLimitsX[2]))
    detectionPositionsY <- c(detectionPositionsY, runif(numFlaseDetections, fieldLimitsY[1], fieldLimitsY[2]))
    detectionIndex <- c(detectionIndex, rep(0, numFlaseDetections))
  }
  detectionList <- data.frame(detectionPositionsX, detectionPositionsY, detectionIndex)
  
  return (detectionList)
}
