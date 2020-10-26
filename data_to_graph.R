data_to_graph <- function(mat, numberOfFrames, arbitrarilySetNoOfPairs) {
  
  vecSta <- seq(
    1, 
    numberOfFrames * arbitrarilySetNoOfPairs, 
    arbitrarilySetNoOfPairs
  )
  vecEnd <- seq(
    arbitrarilySetNoOfPairs, 
    numberOfFrames * arbitrarilySetNoOfPairs, 
    arbitrarilySetNoOfPairs
  )
  diameterVec <- rep(0, numberOfFrames)
  average.path.length.vec.all <- rep(0, numberOfFrames)
  average.path.length.vec.onlyConn <- rep(0, numberOfFrames)
  transitivityVec <- rep(0, numberOfFrames)
  average.degree <- rep(0, numberOfFrames)
  laplacian.spectrum <- list()
  mins <- rep(0, numberOfFrames)
  maxs <- rep(0, numberOfFrames)
  
  for (i in 1:numberOfFrames) {
    
    if (i %% 100 == 0) {
      print(paste(
        "Calculating kernel density of laplacian spectrum for frame:", 
        i
      ))
    }
    
    # ReadFrame
    currFrame <- mat[vecSta[i]:vecEnd[i], ]
    # Create graph
    currGraph <- 
      graph.edgelist(currFrame[which(currFrame[, 1] != 0), ], directed = F)
    # Calculate quantities
    diameterVec[i] <- diameter(currGraph, directed = F)
    average.path.length.vec.all[i] <- 
      average.path.length(currGraph, directed = F, unconnected = F)
    average.path.length.vec.onlyConn[i] <- 
      average.path.length(currGraph, directed = F, unconnected = T)
    average.degree[i] <- mean(degree(currGraph, mode = "total", loops = F))
    transitivityVec[i] <- transitivity(currGraph, type = "undirected")
    # Calculate Laplacian
    adjMat <- get.adjacency(currGraph, type = "both", sparse = F)
    degrees <- rowSums(adjMat)
    laplacian <- diag(nrow(adjMat))
    dummy1 <- which(adjMat == 1, arr.ind = T)[, 1]
    dummy2 <- which(adjMat == 1, arr.ind = T)[, 2]
    laplacian[which(adjMat == 1, arr.ind = T)] <- 
      -(1 / sqrt(degrees[dummy1] * degrees[dummy2]))
    laplacian.spectrum[[i]] <- 
      eigen(laplacian, symmetric = T, only.values = T)[[1]]
    mins[i] <- min(laplacian.spectrum[[i]])
    maxs[i] <- max(laplacian.spectrum[[i]])
  }
  
  return(list(
    diameterVec = diameterVec, 
    average.path.length.vec.all = average.path.length.vec.all, 
    average.path.length.vec.onlyConn = average.path.length.vec.onlyConn, 
    transitivityVec = transitivityVec, 
    average.degree = average.degree, 
    laplacian.spectrum = laplacian.spectrum, 
    mins = mins, 
    maxs = maxs
  ))
}
