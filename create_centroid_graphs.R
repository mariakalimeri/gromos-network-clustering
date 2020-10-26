create_centroid_graphs <- function(mat,
                                   frameNumbersOfCentroids,
                                   numberOfFrames,
                                   arbitrarilySetNoOfPairs,
                                   outname,
                                   clusteringCutoff = clusteringCutoff) {
  vecSta <- 
    seq(1, numberOfFrames * arbitrarilySetNoOfPairs, arbitrarilySetNoOfPairs)
  vecEnd <- 
    seq(
      arbitrarilySetNoOfPairs, 
      numberOfFrames * arbitrarilySetNoOfPairs, 
      arbitrarilySetNoOfPairs
    )
  centroid.graphs <- list()
  # 	for (i in 1:length(frameNumbersOfCentroids)){
  # Write a graph for only the first centroid
  for (i in 1:1) {
    # ReadFrame
    # 		print(i)
    currFrame <- 
      mat[vecSta[frameNumbersOfCentroids[i]]:vecEnd[frameNumbersOfCentroids[i]], ]
    # Create graph
    myG <- graph.edgelist(currFrame[which(currFrame[, 1] != 0), ], directed = F)
    # 		print(paste("Centroid:", i, "done"))
    # ========= Output graph ==============================================
    ## Save graphml
    outname2 <- paste(outname, "centroid", i, sep = "")
    outNameCutoff <- format(round(clusteringCutoff[i], 2), nsmall = 2)
    write.graph(
      myG, 
      file = paste(
        outname2, 
        "_cutoff_", 
        outNameCutoff, 
        ".graphml", 
        sep = ""
      ), 
      format = c("graphml")
    )
    centroid.graphs[[i]] <- myG
  }
  return(centroid.graphs)
}
