# PhyC
Clustering cancer evolutionary trees

## Usage
First, we prepare the list of edge matrices that must include the root and the list of edge length vector.
Main function is phyC. This function have mainly two functionalities. 
  Registration. This includes the function of 
    resolving mono- and multi-furcation tree into bifurcation tree
    completing the number of leaves among the trees, which do not almost affects the calculation of distances
    relabeling considering tree toplogies.

result <- phyC(edgeList,edgeLenList,cluster=3,type='nh',method="ward.D2")

To calculate the diversity of each cluster, we use the diversity function.

result2 <- diversity(resulst)

