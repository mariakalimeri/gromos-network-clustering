calculate_distance_matrix <- function(densitiesObj) {
  distanceMatrix <- matrix(0, length(densitiesObj), length(densitiesObj))
  for (i in 1:(length(densitiesObj) - 1)) {
    if (i %% 100 == 0) {
      print(paste("Difference between", i, "and all the rest"))
    }
    for (j in (i + 1):(length(densitiesObj))) {
      distanceMatrix[i, j] <- 
          sqrt(jensen_shannon(densitiesObj[[i]], densitiesObj[[j]]))
      distanceMatrix[j, i] <- distanceMatrix[i, j]
    }
  }
  return(distanceMatrix)
}
