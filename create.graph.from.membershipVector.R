#rm(list=ls.str())
#library(igraph)

# ====== INPUT ======
## Load membership vectors
### Comment the following line if your matrix is written in a binary file
#mem <- as.matrix(read.table("membershipVector.dat"))
### Comment the following line if your matrix is written in as simple dat file
mem <- as.matrix(read.table("../T260-ala-42/clusteringData/membershipVector_1.00.dat"))

outname <- "networkOfNetworks"
# ===================

# First count the frequencies of each cluster in order to weight the nodes
myHist <- hist(mem, seq(0.5,(max(mem)+0.5), 1), plot=FALSE)
nodeWeights <- myHist$counts
rm(myHist)

# !!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!
# Some older igraph distributions count the vertex indices from 0 and not from 1.
# The lines below are for one-based indices.
# To check that you are fine with your distribution, after creating the myG graph
# object, do min(V(myG))
edgeList <- cbind(mem[1:(length(mem)-1)], mem[2:(length(mem))])
myG <- graph.edgelist(edgeList)

# Simplify the graph by removing first only the loops
myG <- simplify(myG, remove.multiple=FALSE, remove.loops=TRUE)

# Count multiple edges and set them as attribute to the previous graph
edgeWeights <- count.multiple(myG)
newEdgeList <- get.edgelist(myG)
weightedEdgeList <- cbind(newEdgeList, edgeWeights)
simplifiedWeightedEdgeList <- unique(weightedEdgeList)
myG <- simplify(myG, remove.multiple=TRUE)
myG <- set.edge.attribute(myG, "Weight", value=simplifiedWeightedEdgeList[,3])

# Add nodeWeights
V(myG)$Weight <- nodeWeights
myG <- set.vertex.attribute(myG, "Label", value=V(myG))

# # Save the graph as edgeList
#write.table(simplifiedWeightedEdgeList, file=paste(outname,".cg.mcledgelist.txt", sep=""), col.names=F, row.names=F)
## Save graphml with edgeWeights as counts
write.graph(myG, file=paste(outname, ".graphml", sep=""), format=c("graphml"))

## ============================================================================
## ================= NORMALIZATION ============================================
## Build new graph that will have the weights normalized
#simplifiedWeightedEdgeList <- simplifiedWeightedEdgeList[order(simplifiedWeightedEdgeList[,1]),]
#for (i in 1:max(simplifiedWeightedEdgeList)){
#    normConst <- sum(simplifiedWeightedEdgeList[which(simplifiedWeightedEdgeList[,1]==i, arr.ind=T),3])
#    simplifiedWeightedEdgeList[which(simplifiedWeightedEdgeList[,1]==i, arr.ind=T),3] <- simplifiedWeightedEdgeList[which(simplifiedWeightedEdgeList[,1]==i, arr.ind=T),3]/normConst
#}
#myG <- graph.edgelist(simplifiedWeightedEdgeList[,1:2])
#myG <- set.edge.attribute(myG, "Weight", value=simplifiedWeightedEdgeList[,3])
#V(myG)$Weight <- nodeWeights
#myG <- set.vertex.attribute(myG, "Label", value=V(myG))

#write.graph(myG, file=paste(outname, "_normalized.graphml", sep=""), format=c("graphml"))

