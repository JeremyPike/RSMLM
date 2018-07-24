consecutiveFrameGrouping <- function(coords, frames, r) {
  
  numberInGroup = rep(1, dim(coords)[1])

  minFrame = min(frames)

  maxFrame = max(frames) - 1

  for(f in minFrame : maxFrame) {
    
    frameIndices <- which(frames == f | frames == f + 1)
    
    if (length(frameIndices) > 0) {

      coordsTemp <- coords[frameIndices, ]
      framesTemp <- frames[frameIndices]
      numberInGroupTemp <- numberInGroup[frameIndices]

      maxNeigh <- 1
    
      while(maxNeigh > 0) {
        if (!is.vector(coordsTemp)) {
          if(dim(coordsTemp)[1] > 1) {
            fr <- frNN(coordsTemp, r)$id
            numNeigh <- lengths(fr)
            maxNeigh <- max(numNeigh)
            
       
            if (maxNeigh > 0) {
              centreDetection <- which.max(numNeigh)
              groupIndices <- c(centreDetection, fr[[centreDetection]])
              newDetCoords <- colMeans(coordsTemp[groupIndices, ])
              newDetFrame <- f + 1
              newDetNumberInGroup = sum(numberInGroupTemp[groupIndices])
             
              coordsTemp <- coordsTemp[-groupIndices, ]
              coordsTemp <- rbind(coordsTemp, newDetCoords)
              framesTemp <- framesTemp[-groupIndices]
              framesTemp <- c(framesTemp, newDetFrame)
              numberInGroupTemp <- numberInGroupTemp[-groupIndices]
              numberInGroupTemp <- c(numberInGroupTemp, newDetNumberInGroup)
            }
          } else {
            break
          }
          
        } else {
          break
        }
      }
      coords <- coords[-frameIndices, ]
      coords <- rbind(coordsTemp, coords)
      frames <- frames[-frameIndices]
      frames <- c(framesTemp, frames)
      numberInGroup <- numberInGroup[-frameIndices]
      numberInGroup <- c(numberInGroupTemp, numberInGroup)
    }
  }
  detectionList <- cbind(coords, frames, numberInGroup)
  return(detectionList)
}