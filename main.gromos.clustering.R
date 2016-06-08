# This main calls all the other functions in the folder
# 
# WHAT IT DOES:
# Gluster, i.e. group together according to similarity, snapshots of a dynamic 
# network. I wrote and used this script for the output trajectory of a Molecular 
# Dynamics simulation. However, as long as the input format is suitable, the script 
# can be used for any kind of network.
# The input file is a network in the form of an edgelist, i.e. a list of pairs, 
# nodeNumber1 on the left linked to nodeNumber2 on the right. 
# In my case the pairs are oxygens that are hydrogen-bonded to each other in the 
# the hydration shell of a protein. 
# Since there is no fixed number of water molecules across the different trajectory 
# frames, i.e. the hydration shell has variable number of waters since they are moving, 
# the maximum number of pairs is arbitrarily set to some "large" number N (e.g. 3000). 
# Zeros are used to fill it up until line N. So basically the edgelist for 
# frame number 2 starts in line number N+1 and ends in line number 2*N. 
# The edgelist for frame 3 starts in line number 2*N+1 and ends in line number 
# 3*N and so on. The line number N must be given below as an input. See USER INPUT.

# ATTENTION: These scripts will do the calculation for several cases sequentially
# E.g. we are typically interested on the networking of the hydration shell of 
# a certain amino acid at different temperatures.
# For that reason we keep a list of folders starting with a capital "T", 
# one folder for each temperature, and placed in the same level as this code/ folder,
# and start with a (Example name: T380K_histidine). The code enters all folders 
# sequentially, looks for the respective edgelist file named "filename" (see USER INPUT)
# and writes the output in each one of them. 
# !!! Mind to have the same name for all edgelists, each in its own T* folder, 
# and keep in the level of the "code" folder, only the relevant folders starting with
# the letter T (to keep things clean). 

# OUTPUT: The script outputs, inside each T* folder, the membershipVectors 
# (both in dat and plot.pdf forms) and the number of frames that correspond 
# to the centroids (only dat)

rm(list=ls.str())
closeAllConnections()
# ============= START CLOCK ====================================================
ptm <- proc.time()
library(igraph)
# ============= USER INPUT =====================================================
filename <- "pairs.dat"
# ATTENTION: Change the number below to be the arbitrarily set number of pairs
arbitrarilySetNoOfPairs <- 3000
# Choose here cutoff values for clustering. They should be in form of a vector.
# The range of these values can be found by doing some test in a smaller sample 
# of the data.
#clusteringCutoff <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
# Use only one value of clustering on a second run in order to create the networks
# Don't use only one value of the clusteringCutoff if you are not sure whether it 
# create few clusters. One value of clusteringCutoff triggers the creation of centroid
# graphs and it will greate a huge number of files if the cutoff is fine.
clusteringCutoff <- c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
# ============= END INPUT ======================================================
# Loading all the needed scripts. They all must be in the same folder
source(paste(getwd(), "/data.to.graph.R", sep=""))
source(paste(getwd(), "/jensen.shannon.R", sep=""))
source(paste(getwd(), "/bandwidth.calculation.for.density.R", sep=""))
source(paste(getwd(), "/calculate.distanceMatrix.R", sep=""))
source(paste(getwd(), "/gromos.clus.for.probDist.R", sep=""))
source(paste(getwd(), "/create.centroid.graphs.R", sep=""))


# ATTENTION: T* folders must be only the folders that have the input files. 
cases <- system("find ../ -maxdepth 1 -name 'T*'", intern=T)
for (jh in 1:length(cases)){
	print(paste("Inside", cases[jh], "folder"))
	setwd(cases[jh])
	filename <- system(paste("ls -1", filename)), intern=T)
	if (!file.exists("clusteringData")){
    	system("mkdir clusteringData", intern=F) 
    }	
    if (!file.exists("networkAnalysisData")){
    	system("mkdir networkAnalysisData", intern=F) 
    }   
	# ====== DO THE CLUSTERING =====================================================
	# In order to apply the gromos clustering, we need to calculate a distanceMatrix
	# between the density estimates of the laplacian spectrum. 
	# One may want to run the clustering algorithm more than once to try different 
	# cutoffs. Yet the distanceMatrix is needed only once because it is independent
	# of the cutoff. Since its calculation is time consuming, here we first check 
	# if distanceMatrix already exists from a previous run. If yes, then the code
	# moves on to the clustering. If no, then the distanceMatrix will be estimated 
	# again. If you modify parts of the code that are relevant to the distanceMatrix
	# calculation, that is, if you modify any of the following functions
	# -> data.to.graph.R
	# -> jensen.shannon.R
	# -> bandwidth.calculation.for.density.R
	# -> calculate.distanceMatrix.R
	# then you must remove the pre-existing distanceMatrix by removing 
	# the file called distanceMatrix_and_densities.RData.
	if (file.exists("distanceMatrix_and_densities.RData")){
		print("Found existing distanceMatrix. Loading...")
	    load("distanceMatrix_and_densities.RData")
	} else {
	    print("No pre-existing distanceMatrix. Calculating new...")
	    print("Reading input...")
		mat <- as.matrix(read.table(filename))
		numberOfFrames <- nrow(mat)/arbitrarilySetNoOfPairs
	    # Calculate the density estimate of the laplacian spectrum of each frame and 
		# extract it to an object called densitiesObj
		print("Calculating kernel density of laplacian spectrum of graphs...")
		graph.from.mat <- tryCatch(data.to.graph(mat, numberOfFrames=numberOfFrames, arbitrarilySetNoOfPairs=arbitrarilySetNoOfPairs), error=function(e) NULL)
		
    	outname <- substring(cases[jh],4)
    	
    	mean <- mean(graph.from.mat$transitivityVec)
    	std <- sd(graph.from.mat$transitivityVec)
    	write.table(graph.from.mat$transitivityVec, paste("networkAnalysisData/transitivityVec_", outname, ".dat", sep=""), col.names=F, row.names=F)
    
    	mean <- mean(graph.from.mat$diameterVec)
    	std <- sd(graph.from.mat$diameterVec)
    	write.table(graph.from.mat$diameterVec, paste("networkAnalysisData/diameterVec_", outname, ".dat", sep=""), col.names=F, row.names=F)

    	mean <- mean(graph.from.mat$average.path.length.vec.all)
    	std <- sd(graph.from.mat$average.path.length.vec.all)
    	write.table(graph.from.mat$average.path.length.vec.all, paste("networkAnalysisData/average.path.length.vec.all_", outname, ".dat", sep=""), col.names=F, row.names=F)	
    
    	mean <- mean(graph.from.mat$average.path.length.vec.onlyConn)
    	std <- sd(graph.from.mat$average.path.length.vec.onlyConn)
    	write.table(graph.from.mat$average.path.length.vec.onlyConn, paste("networkAnalysisData/average.path.length.vec.onlyConn_", outname, ".dat", sep=""), col.names=F, row.names=F)    
    
    	mean <- mean(graph.from.mat$average.degree)
    	std <- sd(graph.from.mat$average.degree)
    	write.table(graph.from.mat$average.degree, paste("networkAnalysisData/average.degree_", outname, ".dat", sep=""), col.names=F, row.names=F)
		
     	# Use explicitly a different algorithm for the quantile estimation for the bw of the density
		# as the default fails for a large set of our data
		bndwd <- lapply(graph.from.mat$laplacian.spectrum, bandwidth.calculation.for.density) 	
		densitiesObj <- NULL 
		for (i in 1:length(graph.from.mat$laplacian.spectrum)){
		   print(paste("Spectrum:", i))
	 	   densitiesObj[[i]] <- tryCatch(density(graph.from.mat$laplacian.spectrum[[i]], bw=bndwd[[i]])[[2]], error=function(e) NULL)
		}	
		print("Calculating distanceMatrix between laplacian densities...")
	    distanceMatrix <-  tryCatch(calculate.distanceMatrix(densitiesObj), error=function(e) NULL)
	    save(distanceMatrix, densitiesObj, numberOfFrames, file="distanceMatrix_and_densities.RData")
	}	
	# =============== WORKAROUND FOR THE BAD DENSITY ESTIMATIONS ================
	print("Removing problematic frames...")
	dummy <- which(is.nan(distanceMatrix), arr.ind = T)
	if (length(dummy)!=0){
        # NOTE (Tuesday 29/03/2016): When the networks are small (some 10 nodes or less)
        # the density estimates may be wrongly estimated if there exist a few too many
        # eigenvalues equal to 1. The reason is that for the density estimation
        # the quantile calculation is off for too many ones (quantiles collapse 
        # to the same value). Therefore, some densities might be completely off
        # and their respective distances from all other frames will appear as NaN 
        # in the distanceMatrix. 
        # We cure that, simply, by removing the problematic frames. 
        # I need to look only in the first column to find them since they will 
        # repeat in all the columns.
	   	framesToBeRemoved <- dummy[which(dummy[,2]==1),1]
	   	distanceMatrix <- distanceMatrix[-framesToBeRemoved, -framesToBeRemoved]
	   	write.table(framesToBeRemoved, paste("removeFrames.dat"), col.names=F, row.names=F)
	}
	for (i in 1:length(clusteringCutoff)){
		# Do the gromos clustering
		print(paste("Doing the gromos clustering for cutoff", clusteringCutoff[i], "..."))
	    clusGrom <- tryCatch(gromos.clus.for.probDist(distanceMatrix, densitiesObj, cutoff=clusteringCutoff[i]), error=function(e) NULL)
		print("Writing output files...")
		outNameComp <- format(round(clusteringCutoff[i], 2), nsmall = 2)
		write.table(clusGrom$membershipVector, paste("clusteringData/membershipVector_", outNameComp, ".dat", sep=""), col.names=F, row.names=F)
		write.table(clusGrom$centroidPositions, paste("clusteringData/centroidPositions_", outNameComp, ".dat", sep=""), col.names=F, row.names=F)	
		# ====== END DO THE CLUSTERING =================================================
	
		# =========== Get graphs for centroids of leaders ==============================
		## ATTENTION: Use this only when the number of clusters is small, up to 10. It creates and writes down a network in a graphml file 
		## for each centroid. You don't want to do this when you have 200 clusters 
		if (length(clusteringCutoff)==1){
	        print("Reading input connectivity matrix...")
			mat <- as.matrix(read.table(filename))
			print("Writing centroid graphs...")
			outname <- substring(cases[jh],4)
			centroid.graphs <- create.centroid.graphs(mat=mat, frameNumbersOfCentroids=clusGrom$centroidPositions, numberOfFrames=numberOfFrames, arbitrarilySetNoOfPairs=arbitrarilySetNoOfPairs, outname, clusteringCutoff)
		}
		# ======= End get graphs for centroids of leaders ==============================
		
		print(paste("Number of clusters found:",  max(clusGrom$membershipVector)))
		#plot(clusGrom$membershipVector, pch=".", cex=3, col="black", main=paste("Gromos - cutoff=", clusteringCutoff[i], "- # of clus:", max(clusGrom$membershipVector)), xlab="# of frame", ylab="# of 	clusters")
		# =========== PLOT MEMBERSHIPVECTOR=============================================
		print("Plotting membershipVector...")
		if (!is.null(clusGrom)){
			pdf(file=paste("clusteringData/clusterMemberships_", outNameComp, ".pdf", sep=""), onefile=TRUE, paper="special",width=6,height=5)
			#par(mfrow=c(2,3))
			plot(clusGrom$membershipVector, pch=".", cex=3, col="black", xlab="# of frame", ylab="# of clusters", main=paste("Gromos - cutoff=", clusteringCutoff[i], "- # of clus:", max(clusGrom$membershipVector)))
			dev.off()
		}
	}
}
setwd("../code/")
# =========== END CLOCK ========================================================
total.time <- (proc.time() - ptm)/60
print(paste("Total duration in minutes:", total.time[3]))

