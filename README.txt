Glustering, i.e. grouping together according to similarity, snapshots of a dynamic network.

This code was written in the context of a collaboration with Dr. Fabio Sterpone, LBT/IBPC, CNRS, Paris (France)
https://sites.google.com/site/sterponefabio/

Reference: The code was used for the analysis performed in the following publication 
"Configurational Disorder of Water Hydrogen-Bond Network at the Protein Dynamical Transition, Rahaman et all., J Phys Chem B, 2017, http://pubs.acs.org/doi/abs/10.1021/acs.jpcb.7b03888

README file
Author: Maria Kalimeri

=========================================================================================
BASIC INFO AND HOW TO RUN
The main.gromos.clustering.R calls all the other functions in the folder. 
The following description is included in the main.gromos.clustering.R

To use the scripts on your own data file, given that it complies with the 
format requirements (see below), then open main.gromos.clustering.R and 
add/modify, at least, the block between "USER INPUT" and "END INPUT".  

You may use the scripts either in a unix terminal by typing
$ Rscript main.gromos.clustering.R
or from with in R by typing
> source("main.gromos.clustering.R")

===========================================================================================
WHAT DO THE SCRIPTS DO:
Gluster, i.e. group together according to similarity, snapshots of a dynamic 
network. I wrote and used this script for the (modified) output trajectory of a Molecular 
Dynamics simulation. However, as long as the input format is suitable, the script 
can be used for any kind of network.

As the title suggests, clustering is done using the gromos algorithm, described here:
Daura X, Gademann K, Jaun B, Seebach D, van Gunsteren WF, Mark AE. Peptide Folding: 
When Simulation Meets Experiment. Angew Chem Int Ed. 1999;38:236–240. 
doi: 10.1002/(SICI)1521-3773(19990115)38:1/2<236::AID-ANIE236>3.0.CO;2-M. 

While the metric used to compare the distance between the two networks is described here:
"Structural distance and evolutionary relationship of networks"
A. Banerjee, BioSystems, 107 (2012), 186-196
In a nutshell: 
a) Extract the symmetric normalized Laplacian for each network (or snapshot of a dynamic network)
b) Calculate its eigenvalues (i.e. spectrum of the Laplacian)
c) Build a probability density estimate for the spectrum (the spectrum density is said to capture structural characteristics of the network)
d) Use the Jensen-Shannon distance to compare/cluster the probability densities

The input file is a network in the form of an edgelist, i.e. a list of pairs, 
nodeNumber1 on the left linked to nodeNumber2 on the right. 
In my case the pairs are oxygens that are hydrogen-bonded to each other in the 
the hydration shell of a protein. Consequently, some of the description below 
refers to the network of my water molecules. For different cases of networks, just 
read through with an open mind.

(To account for variable node number in input dynamical network)
Since there is no fixed number of water molecules across the different trajectory 
frames, i.e. the hydration shell has variable number of waters since they are moving, 
the maximum number of pairs is arbitrarily set to some "large" number N (e.g. 3000). 
Zeros are used to fill it up until line N. So basically the edgelist for 
frame number 2 starts in line number N+1 and ends in line number 2*N. 
The edgelist for frame 3 starts in line number 2*N+1 and ends in line number 
3*N and so on. The line number N must be given below as an input. See USER INPUT.

ATTENTION: These scripts will do the calculation for several cases sequentially
E.g. we are typically interested on the networking of the hydration shell of 
a certain amino acid at different temperatures.
For that reason we keep a list of folders starting with a capital "T", 
one folder for each temperature, and placed in the same level as this code/ folder,
and start with a (Example name: T380K_histidine). The code enters all folders 
sequentially, looks for the respective edgelist file named "filename" (see USER INPUT)
and writes the output in each one of them. 
!!! Mind to have the same name for all edgelists, each in its own T* folder, 
and keep in the level of the "code" folder, only the relevant folders starting with
the letter T (to keep things clean). 

OUTPUT: The script outputs, inside each T* folder, the membershipVectors 
(both in dat and plot.pdf forms) and the number of frames that correspond 
to the centroids (only dat)
