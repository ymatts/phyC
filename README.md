# PhyC
R packages for clustering cancer evolutionary trees. 

Depend:ape, igraph

## Usage
First, we prepare the list of edge matrices that must include the root and the list of edge length vector. The edge matrix represents cancer evolutionary topology and the edge length represents the number of additional SSNVs from the parental suc-clone. Here is illustration.

<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/ssnv.png" width="600" height="400" />


Main function is phyC. This function have mainly two functionalities. 
* Registration. 
    + This includes the function of 
        + resolving mono- and multi-furcation tree into bifurcation tree
        + completing the number of leaves among the trees, which do not almost affects the calculation of distances
        + relabeling considering tree toplogies.
    - The overall scheme is illustrated here.
     <img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/regist.png" width="600" height="400" />

    - Here are some example of registrations. The example shows the mono-furcation and multi-furcation trees (upper left most and lower left most). The resoved toplogies are shown in the next two. The upper & lower right most show the illustration of completing the number of leaves (in this case, 64 leaves). This procedure is need to deal with trees in tree space (Birella,et al. 2001).   
<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/regis_example.jpeg" width="600" height="400" />

* Clustering trees.
    + Hierarchical(Ward's method)
    + Non-hierarchical(k-means)
<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/clust.png" width="600" height="400" />


We also provide the function for interpretating resulting clusters. The diversity is defined asthe number of sub-clones at a SSNVs accumulations. This gives us the insights for accerelation of sub-clonal expansions. 
<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/diversity.png" width="600" height="400" />

######Use of phyC
The phyC needs the edgeList, edgeLenList and cluster(the number of the cluster) in minimal. Here is an example.

```r:phyC.R
result <- phyC(edgeList,edgeLenList,cluster=3)
```

######Use of diversity
To calculate the diversity of each cluster, we use the diversity function. This function requires only the phyC object.

```r:diversity.R
result2 <- diversity(result)
```

##Example
