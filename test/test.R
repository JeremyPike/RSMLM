#!/usr/bin/env Rscript

library(RSMLM)
library(pryr)

#### check algorithm is reproducible ####

# clustering parameters
r <- 20
thresh <- 10

# load example data
detectionList <- read.csv('detectionList.csv')
# extract x and y coordinates
coords <- as.matrix(detectionList[, c('x', 'y')])
# perform clustering
clusterIndices <- clusterTomato(coords, r, thresh)

# check it produces the same results
check <- all(clusterIndices == detectionList$tom)

if (check) {
  print('Reproduces result')
} else {
  print('Does not reproduce result')
}


#### Demonstrate memory leak #####

numRepeats <- 1000
numPoints <- 10000

# generate some randomly placed 2D coordinates
coords <- matrix(0, numPoints, 2)
for (d in 1 : 2) {
  coords[, d] <- runif(numPoints, c(0, 2000), c(0, 2000))
}

# repeat clustering to domstrate memory build up
mem_base <- mem_used()
mem_prev <- 0.0
# vector = c() # this is just to check if the memory print works as intended
cat(sprintf("Initial memory %s Mb \n", mem_base*1e-6))
for (n in 1 : numRepeats) {
  
  clusterIndices <- clusterTomato(coords, 30, 10) # replaced with the expanded code
  
  # ############################
  # # Expanded code

  # library("dbscan")

  # r <- 30
  # threshold <- 10
  # coords2 <- as.matrix(coords)
  
  # numDimensions <- dim(coords2)[2]
  
  # # compute density using number of other detections within fixed search radius (library dbscan)
  # density <- pointdensity(coords2, r, type = "frequency")
  # # tomato algorithm performed using c++ function
  # clusterIndices <- tomatoDens(coords2, density, r, threshold)$clusters # <<-- here's the problem !

  # # Expanded code --END
  # ############################

  mem_current <- mem_used() - mem_base
  mem_diff <- mem_current - mem_prev
  cat(sprintf("Iteration %i additional memory used: %f Mb (increased of %f Mb from the previous iteration)\n", n, mem_current*1e-6, mem_diff*1e-6))
  mem_prev <- mem_current
}
  