# phyC: Clustering Cancer Evolutionary Trees
The R package for clustering cancer evolutionary trees. 

<strong>Version:</strong>
0.2.0 (2016.08.01)

<strong>Depends:</strong>

R(>=3.2.2)

ape, igraph, ggplot2, grDevices, png, RColorBrewer

<strong>Authour:</strong>

Yusuke Matsui

<strong>Contact:</strong>

ymatsui[at]med.nagoya-u.ac.jp


## General overview
Multi-regional sequencing provide new opportunities to investigate genetic heterogeneity within- or between- common tumours from a view of evolutionary perspective. The state-of-the-art methods have been proposed for reconstructing cancer sub-clonal evolutionary trees based on multi-regional sequencing data to develop models of cancer evolution. However, the methods developed thus far are not sufficient to characterize and interpret the diversity of cancer sub-clonal evolutionary trees. We propose a clustering method (phyC) for cancer sub-clonal evolutionary trees, in which sub-groups of the trees are identified based on topology and edge length attributes. For interpretation, we also propose a method for evaluating the diversity of trees in the clusters, which provides insight into the acceleration of sub-clonal expansion. 

## Introduction

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

Main outpus are
* cluster assignments
* Multidimensional scaling configuration of trees in clusters
* sub-clonal diversity of trees in clusters

In case of VAF input, we need reconstruct cancer sub-clonal evolutionary trees using existing methods. We implement two method; maximum parsimony approch and clustering-based rooted-constraint network approach. The former is based on acctran in phangorn package (Klein, et al. 2010) and the latter is based LICHeE (Popic, et al. 2015). 

The meaning of each element of reconstructed trees is as follows: the root and its subsequent node represent a normal cell and founder cell, respectively. Sub-clones are described as nodes below the founder cell, and edge lengths indicate the number of SSNVs that are newly accumulated in descendant nodes.


## Package overview

####phyC
Main function of this package is <em>phyC</em>. This function has mainly the two functionalities; registration and clustering. 

######Registration
* This includes; 
 + resolving mono- and multi-furcation tree into bifurcation tree
 + completing the number of leaves among the trees, which do not almost affects the calculation of distances
 + relabeling to remove the label differences so that the same trees with the same toplogies are regarded as the same labelled trees.
 + normalizing edge lengths

* The reason why we need the registration above are listed point by point;
 + The <em>tree space</em> assumes the <em>n-tree</em> (strictly bifurcation) but actual data do not take such forms
 + The number of leaves must be the same in the <em>tree space</em>, but actually not.
 + The labels of leaves are distinguishied in the <em>tree space</em>, but the patterns of the accumulated SSNVs are rarely identical among the patients since there are too many patterns of SSNVs accumulations in the actual cases (we cannot divide the patients into several sub groups with the labels). Thus we focus on the patters of the number of accumulated SSNVs, that is, we regard the trees with the same toplogies as the same labelled trees. 
 + The number of accumulated SSNVs are also different from patients to patients, or studies to studies because of sequencing depth, errors and so on. 

The overall scheme of the registration is illustrated here (Figure 3). We at first prepare the maximal trees and then we encode the observed tree toplogies(from root to leaf, left to right). The collapsed edges are regarded as zero length edges. The solid lines and dotted lines indicate the encoded toplogies and the collapsed edges, respectively.
    
<center> <img src="https://github.com/ymatts/PhyC/blob/master/img/regist.png" alt="overall scheme of the registration" width="600" height="400"><center>

<strong>Figure 3.</strong> Overall scheme of the registration

Here are some example of registrations (Figure 4). The example shows the mono-furcation and multi-furcation trees (upper and lower left panels). The resoved toplogies are shown in the next two sides (2nd and 3rd panels from the left). The upper and lower right most panels show the illustration of completing the number of leaves (in this case, 64 leaves). PhyC automatically performs these procedures.

<img src="https://github.com/ymatts/PhyC/blob/master/img/regis_example.jpeg" align="center" alt="an example of the registration"  width="600" height="400">

<strong>Figure 4.</strong> Example of the registration

The second functionality of <em>phyC</em> is the clustering.

######Clustering trees
   + Hierarchical (Ward's method)
   + Non-hierarchical (k-means)

Those clustering algorithms are naturally extended from Euclidean space to tree space using the tree distances. Here is an example of non-hierarchical clustering with 4 clusters (Figure 5).

<img src="https://github.com/ymatts/PhyC/blob/master/img/phyC.plot.png" align="center" width="600" height="400" />

<strong>Figure 5.</strong> Example of non-hierarchical clustering

####diversity
The next function <em>diversity</em> is for interpretating clusters. We need some quantifications of trees in the clusters to evaluate them efficiently. We define the diversity as the number of sub-clones at a levels of SSNVs accumulations. In Figure 6, <em>h</em> indicates the normalized number of additional SSNVs (ranged [0,1]) and <em>g(h)</em> indicates the the number of subclones. This gives us the insights for accerelation of sub-clonal expansions. 

<img src="https://github.com/ymatts/PhyC/blob/master/img/diversity.plot2.png" align="center" width="600" height="400" />

<strong>Figure 6.</strong> Example of diversities for clusters

####phyCMD
We also provide the function <em>phyCMD</em> for plotting configurations of trees as shown in Figure 2.

####lichee2edge
As a utility, we implement the function <em>lichee2edge</em> with which we can obtain the evloutionary trees from variant allele frequencies by LICHeE (Popic,et al. 2015).

##Usage
######Installation

```r:install_git.R
devtools::install_git(url = "https://github.com/ymatts/PhyC")
```

######Use of phyC
The phyC needs the edgeList, edgeLenList and cluster(the number of the cluster) in minimal. 

```r:phyC.R
result <- phyC(edgeList,edgeLenList,cluster=3)
```

Here is an example.
```r:phyC.R
library(PhyC)
data(evol)
result <- phyC(evol$edgeList,evol$edgeLenList,cluster=4,type='nh')
result$cluster # Assignment of clusters
phyC.plot(result) # Output the plot as Figure 5.
```


######Use of diversity
To calculate the diversity of each cluster, we use the diversity function. This function requires only the phyC object.

```r:diversity.R
result2 <- diversity(phyC.obj)
```

Here is an example.

```r:diversity.R
library(PhyC)
data(evol)
result <- phyC(evol$edgeList,evol$edgeLenList,cluster=4,type='nh')
result2 <- diversity(result) # Output the plot as Figure 6.
```

######Use of phyCMD
To obtain the configuration of trees in Euclidean space as Figure 2, we use phyCMD function. The input is phyC object.

```r:phyCMD.R
result3 <- phyCMD(phyC.obj)
```


```r:phyCMD.R
library(PhyC)
data(evol)
result <- phyC(evol$edgeList,evol$edgeLenList,cluster=4,type='nh')
result3 <- phyCMD(result) # Output the plot as Figure 6.
```


######Use of lichee2edge
To reconstruct the cacer evolutionary trees, we adopt LICHeE (Popic, et al. 2015). You need to get the LICHeE engine from <a href="url">http://viq854.github.io/lichee/</a>. Here we provide the utility function to utilize the LICHeE from R. The input is VAF matrix (the format is described in the above site). You can specify the parameters of LICHeE as you need. If you don't specify them, lichee2edge automatically set the default values that are suggested in (Popic et al. 2015). The output is edge matrix and edge length. 


```r:lichee2edge.R
result4 <- lichee2edge('Path to lichee.jar',vaf)
```

```r:phyCMD.R
library(PhyC)
data(vaf) #The list of VAF matrices of ccRCC and HGSC studies

#### We set the parameter of absent, present, minPrivateClsuterSize, maxClusterDist for each patient.

absent <- c(rep(0.005,6),0.01,0.005,rep(0.005,4),0.01,0.005)
present <- c(rep(0.005,6),0.01,0.005,rep(0.01,4),0.04,0.01)
minPrivateClusterSize  <- c(rep(1,7),2,rep(1,3),3,3,1,2)
maxClusterDist <- c(rep(0.2,9),rep(0.2,5),0.1)

result4 <- vector("list",length(vaf))
for(i in seq_along(vaf)){  
  licheeParamIO <- list(normal=1) ## INPUT/OUTPUT AND DISPLAY OPTIONS
  licheeParamFilter <- list(absent=absent[i],present=present[i]) ## SSNV FILTERING AND CALLING parameter
  licheeParamPhy <- list(minPrivateClusterSize=minPrivateClusterSize[i],maxClusterDist=maxClusterDist[i]) ## PHYLOGENETIC NETWORK CONSTRUCTION AND TREE SEARCH parameter
  result4[[i]] <- lichee2edge("LICHeE/release",vaf[[i]],licheeParamIO,licheeParamFilter,licheeParamPhy)
}
```

## References
1. Matsui Y, Niida A, Uchi R. Mimori K, Miyano S, and Shimamura T.(2016) Clustering cancer evolutionary trees. (submitted). 
2. Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S.(2015) Fast and scalable inference of multi-sample cancer lineages. Genome Biol. 16:91.
3. Beerenwinkel N, Schwarz RF, Gerstung M, Markowetz F. (2014)  Cancer  evolution:   mathematical  models  and computational inference. Syst Biol. 6(1):e-2
4. Billera LJ, Holmes SP, Vogtmann K. (2001) Geometry of the Space of Phylogenetic Trees. Adv. Appl.Math. 27(4),733-767.
