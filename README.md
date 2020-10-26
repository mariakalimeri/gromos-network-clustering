# Clustering snapshots of a dynamic network

A collection of scripts for clustering, i.e. grouping together according to 
similarity, snapshots of a dynamic network.

This code was written in the context of a collaboration with [Dr. Fabio Sterpone, 
LBT/IBPC, CNRS, Paris (France)](https://sites.google.com/site/sterponefabio/).

The code was used for the analysis performed in the following publication 
["Configurational Disorder of Water Hydrogen-Bond Network at the Protein 
Dynamical Transition, Rahaman et all., J Phys Chem B, 2017](http://pubs.acs.org/doi/abs/10.1021/acs.jpcb.7b03888)

## Usage 

The `main_gromos_clustering.R` script calls all the other functions in the 
folder. 

To use the scripts on your own data file, given that it complies with the 
format requirements (see below), open `main_gromos_clustering.R` and 
add/modify, at least, the block between `USER INPUT` and `END INPUT`.  

Run the script from within R

```r
> source("main_gromos_clustering.R")
```

or from a unix shell

```
Rscript main_gromos_clustering.R
```

### Details 

The scipts cluster together, according to similarity, snapshots of a dynamic 
network. I wrote and used this script for the (modified) output trajectory of a 
Molecular Dynamics simulation. However, as long as the input format is suitable, 
the script can be used for any kind of network.

Clustering is done using the gromos algorithm, as described in
["Peptide Folding: When Simulation Meets Experiment", Daura X et al. 
Angew Chem Int Ed.
1999;38:236–240](https://onlinelibrary.wiley.com/doi/10.1002/%28SICI%291521-3773%2819990115%2938%3A1/2%3C236%3A%3AAID-ANIE236%3E3.0.CO%3B2-M).

The metric used to compare the distance between the two networks is described 
here: ["Structural distance and evolutionary relationship of networks"
A. Banerjee, BioSystems, 107 (2012), 186-196](https://www.sciencedirect.com/science/article/abs/pii/S0303264711001869).

In a nutshell: 

1. Extract the symmetric normalized Laplacian for each network (or snapshot of 
a dynamic network)
2. Calculate its eigenvalues (i.e. spectrum of the Laplacian)
3. Build a probability density estimate for the spectrum (the spectrum density 
is said to capture structural characteristics of the network)
4. Use the Jensen-Shannon distance to compare/cluster the probability densities

The input file is a network in the form of an edgelist, i.e. a list of pairs, 
node_1 on the left linked to node_2 on the right. 
In my case the pairs are oxygens that are hydrogen-bonded to each other in the 
the hydration shell of a protein. Consequently, some of the description below 
refers to the network of my water molecules. For different cases of networks, just 
read through with an open mind.

To account for the variable number of nodes in the input network,
since there is no fixed number of water molecules across the different 
trajectory frames -- i.e. in the hydration shell some waters move in while some 
move out -- the maximum number of pairs is arbitrarily set to some "large" 
number `N` (e.g. 3000). Zeros are used to fill it up until line `N`. That is,
the edgelist for frame number 2 starts in line number `N+1` and ends 
in line number `2*N`. The edgelist for frame 3 starts in line number `2*N+1` and 
ends in line number `3*N` and so on. The line number `N` must be given below as 
an input. See `USER INPUT`.

> Important note: In their current version, the scripts will do the calculation 
for several cases sequentially. E.g. for the water network clustering case, 
we are interested on the networking of the hydration shell, around a certain 
amino acid, at different temperatures. For that purpose, we need a list of 
subfolders, one for each temperature, starting with a capital `T*` and placed 
in the root directory of this repository. An example name is `T380K_histidine` 
(The `T` prefix of the folder name is actually hard-coded. You'd need to modify 
the code below the `USER INPUT` if you want to change that). The code will 
parse all folders with the `T*` prefix and look for the respective edgelist
file, with the name `filename` as defined in the `USER INPUT` section. 
The respective output is written in each 
subfolder. Keep the same `filename` for the edgelists in each `T*` subfolder. 

### Code output 

The script outputs, inside each `T*` subfolder, the membershipVectors 
(both in a `.dat` format as well as a `plot.pdf` file) and the number of frames 
that correspond to the centroids (this is only a `.dat` file)
