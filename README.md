# PhyC
R packages for clustering cancer evolutionary trees

Depend:ape, igraph

## Usage
First, we prepare the list of edge matrices that must include the root and the list of edge length vector.
Main function is phyC. This function have mainly two functionalities. 
* Registration. 
    + This includes the function of 
    + resolving mono- and multi-furcation tree into bifurcation tree
    + completing the number of leaves among the trees, which do not almost affects the calculation of distances
    + relabeling considering tree toplogies.
* Evaluating diversity. 
    + This function evaluates the accerelation of sub-clonal expansion of trees. Graphical ouput is also supported.

######Use of phyC
The phyC need the edgeList, edgeLenList and cluster(the number of the cluster) in minimal. Here is an example.

```r:phyC.R
result <- phyC(edgeList,edgeLenList,cluster=3)
```

######Use of diversity
To calculate the diversity of each cluster, we use the diversity function. This function requires only the phyC object.

```r:diversity.R
result2 <- diversity(result)
```
