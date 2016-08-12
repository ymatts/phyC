# phyC: Clustering Cancer Evolutionary Trees
The R package for clustering cancer evolutionary trees. 

<strong>Version:</strong>
0.2.0 (2016.08.01)

<strong>Depends:</strong>

R(>=3.2.2)

ape, igraph, ggplot2, grDevices, png, RColorBrewer

<strong>Author:</strong>

Yusuke Matsui

<strong>Contact:</strong>

ymatsui[at]med.nagoya-u.ac.jp


## General overview
Multi-regional sequencing provide new opportunities to investigate genetic heterogeneity within- or between- common tumours from a view of evolutionary perspective. The state-of-the-art methods have been proposed for reconstructing cancer sub-clonal evolutionary trees based on multi-regional sequencing data to develop models of cancer evolution. However, the methods developed thus far are not sufficient to characterize and interpret the diversity of cancer sub-clonal evolutionary trees. We propose a clustering method (phyC) for cancer sub-clonal evolutionary trees, in which sub-groups of the trees are identified based on topology and edge length attributes. For interpretation, we also propose a method for evaluating the diversity of trees in the clusters, which provides insight into the acceleration of sub-clonal expansion. 

## Input/Output and breif description of phyC

phyC (Phylogenetic tree Clustering) is designed for classifying cancer evolutionary trees. 

The main inputs of phyC is
* Variant allele frequency (VAF)

or alternatively

* Edge matrix and edge length vector from reconstructed sub-clonal evolutionary trees

<strong>Table 1.</strong> Example of VAF
<table>
<tr> <th>  </th> <th> Normal </th> <th> Region 1 </th> <th> Region 2 </th> <th> Region 3 </th> <th> Region 4 </th> <th> Region 5 </th>  </tr>
  <tr> <td align="right">gene 1 </td> <td align="right"> 0.00 </td> <td align="right"> 0.24 </td> <td align="right"> 0.18 </td> <td align="right"> 0.24 </td> <td align="right"> 0.24 </td> <td align="right"> 0.22 </td> </tr>
  <tr> <td align="right">gene 2 </td> <td align="right"> 0.00 </td> <td align="right"> 0.15 </td> <td align="right"> 0.08 </td> <td align="right"> 0.12 </td> <td align="right"> 0.25 </td> <td align="right"> 0.24 </td> </tr>
  <tr> <td align="right">gene 3 </td> <td align="right"> 0.00 </td> <td align="right"> 0.19 </td> <td align="right"> 0.18 </td> <td align="right"> 0.22 </td> <td align="right"> 0.26 </td> <td align="right"> 0.26 </td> </tr>
   </table>
   
   
<strong>Table 2.</strong> Example of edge matrix and corresponding edge length vector
<table>
<tr> <th>  </th> <th>Node 1 </th> <th>Node 2 </th>  </tr>
  <tr> <td align="right">Edge 1 </td> <td align="right"> 9 </td> <td align="right"> 8 </td> </tr>
  <tr> <td align="right">Edge 2 </td> <td align="right"> 8 </td> <td align="right"> 6 </td> </tr>
  <tr> <td align="right">Edge 3 </td> <td align="right"> 6 </td> <td align="right"> 7 </td> </tr>
  <tr> <td align="right">Edge 4 </td> <td align="right"> 6 </td> <td align="right"> 4 </td> </tr>
  <tr> <td align="right">Edge 5 </td> <td align="right"> 4 </td> <td align="right"> 1 </td> </tr>
  <tr> <td align="right">Edge 6 </td> <td align="right"> 1 </td> <td align="right"> 3 </td> </tr>
  <tr> <td align="right">Edge 7 </td> <td align="right"> 8 </td> <td align="right"> 10 </td> </tr>
  <tr> <td align="right">Edge 8 </td> <td align="right"> 4 </td> <td align="right"> 5 </td> </tr>
  <tr> <td align="right">Edge 9 </td> <td align="right"> 1 </td> <td align="right"> 2 </td> </tr>
   </table>

<table>
<tr> <th>Edge 1 </th> <th>Edge 2 </th> <th>Edge 3 </th> <th>Edge 4 </th> <th>Edge 5 </th> <th>Edge 6 </th> <th>Edge 7 </th> <th>Edge 8 </th> <th>Edge 9 </th>  </tr>
  <tr> <td align="right"> 209.00 </td> <td align="right"> 0.00 </td> <td align="right"> 17.00 </td> <td align="right"> 0.00 </td> <td align="right"> 8.00 </td> <td align="right"> 2.00 </td> <td align="right"> 44.00 </td> <td align="right"> 44.00 </td> <td align="right"> 19.00 </td> </tr>
   </table>

In case of VAF input, we need reconstruct cancer sub-clonal evolutionary trees using existing methods. We implement two method; maximum parsimony approch and clustering-based rooted-constraint network approach. The former is based on acctran in phangorn package (Klein, et al. 2010) and the latter is based LICHeE (Popic, et al. 2015). 

The meaning of each element of reconstructed trees is as follows: the root and its subsequent node represent a normal cell and founder cell, respectively. Sub-clones are described as nodes below the founder cell, and edge lengths indicate the number of SSNVs that are newly accumulated in descendant nodes.

phyC peform clustering for a set of the reconstructed evolutionary trees after transforming tree objects via transformation of the tree topologies and edge attributes to allow for effective comparison among trees, a procedure we refer to as tree registration.

Main outpus are
* Cluster assignments
* Multidimensional scaling (MDS) configuration of trees in clusters
* Sub-clonal diversity of trees in clusters

Clustering is based on standard hierarchical clustering algorithm and default algorithm of phyC is Ward's method. 
For visualization of clustering results for interpretation, we develop the two tools; MDS and sub-clonal diversity plot. In MDS, we approximately embed the registered trees into lower-dimensional Euclidean space. We overlay tree shapes over embedded Euclidean space in the plot. In sub-clonal diversity plot, we apply the concept of a lineage-through-time (LTT) plot to visualize how sub-clones evolve along SSNV accumulation. The LTT plot generally describes the time vs. number of lineages; and in the present case, this is expressed as the number of sub-clones (y-axis) vs. the fraction of accumulated SSNVs (x-axis).


##Usage
#####Installation

```r:install_git.R
library(devtools)
install_git(url = "https://github.com/ymatts/PhyC")
```

#####Tree reconstruction from VAF profile
The phyC implement two reconstruction methods from VAF profiles: Maximum parsimony method implemented as acctran in phangorn package by (Klaus, 2011) and LICHeE (Popic et al. 2015).

######par.tree
```r:par_tree.R
tree <- par.tree(VAF, thr = 0.05) ## maximum parsimony method. 
```

An example is here.

```r:par_tree.R
library(phyC)
data(ccRCC)
vaf <- lapply(ccRCC,function(x)x[,-(1:3)]) ## par.tree only needs a matrix of gene * (Normal / Region VAF).
label <- names(vaf)
trees <- vector("list",length(vaf))
for(i in seq_along(vaf)){
  trees[[i]] <- par.tree(vaf[[i]])
}
edgeList <- lapply(trees,function(x)x$edge) ## obtain edge list
edgeLenList <- lapply(trees,function(x)x$edge.length) ## obtain edge length list
```

######lichee2edge
The other method is LICHeE (Popic et al. 2015). You need to get the LICHeE engine from <a href="url">http://viq854.github.io/lichee/</a>. Here we provide the utility function to utilize the LICHeE from R. The input is VAF matrix (the format is described in the above site). You can specify the parameters of LICHeE as you need. If you don't specify them, lichee2edge automatically set the default values that are suggested in (Popic et al. 2015). The output is edge matrix and edge length. 
```r:lichee2edge.R
tree <- lichee2edge('Path to lichee.jar', VAF, parameters)
```
```r:lichee2edge.R
absent <- c(rep(0.005,6),0.01,0.005)
present <- c(rep(0.005,6),0.01,0.005)
minPrivateClusterSize  <- c(rep(1,7),2)
maxClusterDist <- c(rep(0.2,8))

trees <- vector("list",length(vaf))
for(i in seq_along(vaf)){  
  licheeParamIO <- list(normal=1) ## INPUT/OUTPUT AND DISPLAY OPTIONS
  licheeParamFilter <- list(absent=absent[i],present=present[i]) ## SSNV FILTERING AND CALLING parameter
  licheeParamPhy <- list(minPrivateClusterSize=minPrivateClusterSize[i],maxClusterDist=maxClusterDist[i]) ## PHYLOGENETIC NETWORK CONSTRUCTION AND TREE SEARCH parameter
  trees[[i]] <- lichee2edge("LICHeE/release",ccRCC[[i]],licheeParamIO,licheeParamFilter,licheeParamPhy)
  edgeList <- lapply(trees,function(x)x$edge) ## obtain edge list
  edgeLenList <- lapply(trees,function(x)x$edge.length) ## obtain edge length list
}
```
#####Registrating trees and Clustering
######phyC
The main function phyC needs the edgeList, edgeLenList and cluster(the number of the cluster) in minimal. 

```r:phyC.R
result <- phyC(edgeList,edgeLenList,cluster)
```

Here is an example. We use the edgeList and edgeLenList from maximum parsimony method here.
```r:phyC.R
result <- phyC(edgeList,edgeLenList,cluster=3,type='h')
result$cluster # Assignment of clusters
```


######sub-clonal diversity plot
To calculate the diversity of each cluster, we use the diversity function. This function requires only the phyC object.
```r:diversity.R
result2 <- diversity(phyC.obj)
```

Here is an example.

```r:diversity.R
result2 <- diversity(result)
```

######phyCMD
To obtain the configuration of trees in Euclidean space, we use phyCMD function. The input is phyC object.

```r:phyCMD.R
result3 <- phyCMD(phyC.obj)
```

Here is an example.

```r:phyCMD.R
result3 <- phyCMD(result) # Output the plot as Figure 6.
```

## References
1. Beerenwinkel N, Schwarz RF, Gerstung M, Markowetz F. (2014)  Cancer  evolution:   mathematical  models  and computational inference. Syst Biol. 6(1):e-2
2. Billera LJ, Holmes SP, Vogtmann K. (2001) Geometry of the Space of Phylogenetic Trees. Adv. Appl.Math. 27(4),733-767.
3. Klaus,P.S. (2011) phangorn: phylogenetic analysis in R. Bioinformatics, 27, 592-593.
4. Matsui Y, Niida A, Uchi R. Mimori K, Miyano S, and Shimamura T.(2016) Clustering cancer evolutionary trees. (submitted, preprint version: BIORXIV/2016/069302).
5. Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S.(2015) Fast and scalable inference of multi-sample cancer lineages. Genome Biol. 16:91.
