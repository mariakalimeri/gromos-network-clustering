gromos.clus.for.probDist <- function(distanceMatrix, densitiesObj, cutoff){
    
    centroid <- list()
    centroidVector <- NULL
    membershipVector <- rep(0, numberOfFrames)  
    indexCluster <- 0
    yesNoMat <- (distanceMatrix <= cutoff)
    rowSummation <- rowSums(yesNoMat)
    numberOfFrames <- length(rowSummation)
    clusteredConformationsSoFar <- 0
    while (max(rowSummation)!=0){
        clusteredConformationsSoFar <- clusteredConformationsSoFar + max(rowSummation)
        print(paste(100*(clusteredConformationsSoFar/numberOfFrames), "% of conformations clustered"))
        indexCluster <- indexCluster+1
        whichGuy <- which(rowSummation==max(rowSummation))[1]
        toBeZero <- c(whichGuy, which(yesNoMat[whichGuy,]==1))
	    centroid[[indexCluster]] <- densitiesObj[[whichGuy]]
        centroidVector <- c(centroidVector, whichGuy)
        membershipVector[toBeZero] <- indexCluster 
        yesNoMat[toBeZero,] <- 0
        yesNoMat[,toBeZero] <- 0
        rowSummation <- rowSums(yesNoMat)
    }
    return(list(numberOfClusters=max(membershipVector), membershipVector=membershipVector, centroidDistribution=centroid, centroidPositions=centroidVector))
}
