# PhyC
R packages for clustering cancer evolutionary trees. 

Depend:ape, igraph

## General overview

Recently a lot of the reconstruction methods of cancer sub-clonal evolutionary trees based on the variant allele frequency are proposed. Those methods can automatically produce the evolutionary trees. However, interpretation of the large number of trees from those methods remain unresolved. We propose the classification method for the cancer sub-clonal evolutionary trees. In this method, we adopt the tree space theory to analyse the tree set. The proposed method can identify the cluster structure based on the tree topology and edge length attributes. For the interpretation of the clusters, we also provide the method for evaluating the diversity of trees in the clusters, which gives an insight for the acceleration of sub-clonal expansion.

## Introduction

PhyC (Phylogenetic tree Clustering) is designed for classifying cancer evolutionary trees. 

The main inputs are
* list of edge matrices that must include the root
* list of edge length vector

and the main outpus are
* cluster assignments
* registered trees
* accerelation of sub-clonal expansions

In the input, the edge matrix represents cancer evolutionary topology and the edge length represents the number of additional Somatic Single Nucleotide Variants (SSNVs) from the parental suc-clone. Here is an illustration (Figure 1).

<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/ssnv.png" width="600" height="400" />

<strong>Figure 1.</strong>

Our model assumes cancer clonal theory (Nowell,1976; Nik-Zainal {\it et al.,} 2012) that
 
1. no mutation occurs twice in the course of cancer evolution
2. no mutation is ever lost.

Our classification method takes the inputs from those reconstruction methods below that outputs the various cancer evolutionary trees from the multiregional cancer sequencing per patiant and the excellent overview of those methods are described in (Beerenwinkel {\it et al.,} 2014); PhyloSub (Jiao {\it et al.,} 2013), PyClone (Roth {\it et al.,} 2014), SciClone (Miller {\it et al.,} 2014), Clomial (Zare {\it et al.,} 2014), Trap (Strino {\it et al.,} 2013), SubcloneSeeker (Qiao {\it et al.,} 2014), and LICHeE (Popic {\it et al.,} 2015) etc.

PhyC classifies those evolutinary trees based on the tree toplogies and edge length attributes under the framework of <em>tree space</em> (Birella,et al. 2001). Here is an example of a result (Figure 2).

<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/mds.png" width="600" height="400" />

<strong>Figure 2.</strong> Example of the clustering

## Package overview

####phyC
Main function of this package is <em>phyC</em>. This function have mainly the two functionalities; registration and clustering. 

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

The The overall scheme of the registration is illustrated here (Figure 3). We at first prepare the maximal trees and then we encode the observed tree toplogies(from root to leaf, left to right). The collapsed edges are regarded as zero length edges. The solid lines and dotted lines indicate the encoded toplogies and the collapsed edges, respectively.
    
<center> <img src="https://github.com/ymatts/PhyC/blob/master/img/regist.png" alt="overall scheme of the registration" width="600" height="400"><center>

<strong>Figure 3.</strong> Overall scheme of the registration

Here are some example of registrations (Figure 4). The example shows the mono-furcation and multi-furcation trees (upper and lower left panels). The resoved toplogies are shown in the next two sides (2nd and 3rd panels from the left). The upper and lower right most panels show the illustration of completing the number of leaves (in this case, 64 leaves). PhyC automatically performs these procedures.

<img src="https://github.com/ymatts/PhyC/blob/master/img/regis_example.jpeg" align="center" alt="an example of the registration"  width="600" height="400">

<strong>Figure 4.</strong> Example of the registration"

The second functionality of <em>phyC</em> is the clustering.

######Clustering trees
   + Hierarchical(Ward's method)
   + Non-hierarchical(k-means)

Those clustering algorithms are naturally extended from Euclidean space to tree space using the tree distances. Here is an example of non-hierarchical clustering with 4 clusters (Figure 5).

<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/clust.png" width="600" height="400" />

<strong>Figure 5.</strong> Example of non-hierarchical clustering

####diversity
The secondary function <em>diversity</em> is for interpretating clusters. We need some quantifications of trees in the clusters to evaluate them efficiently. We define the diversity as the number of sub-clones at a levels of SSNVs accumulations. This gives us the insights for accerelation of sub-clonal expansions. 
<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/diversity.png" width="600" height="400" />

<strong>Figure 6.</strong> Exampmle of diversities for clusters

####phyMDS
We also provide the function <em>phyMDS</em> for plotting configurations of trees as shown in Figure 2. This simply performed with classical multidimensional scaling.

####lichee2edge
As a utility, we implement the function <em>lichee2edge</em> with which we can obtain the evloutionary trees from variant allele frequencies by LICHeE (Popic,2015).

##Usuage

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
