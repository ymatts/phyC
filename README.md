# PhyC: Clustering Cancer Evolutionary Trees
The R package for clustering cancer evolutionary trees. 

<strong>Version:</strong>
0.1.2 (2016.1.18)

<strong>Depends:</strong>

R(>=3.2.2)

ape, igraph, ggplot2, grDevices, png, RColorBrewer

<strong>Authour:</strong>

Yusuke Matsui

<strong>Contact:</strong>

ymatsui[at]med.nagoya-u.ac.jp


## General overview

Recently a lot of the reconstruction methods of cancer sub-clonal evolutionary trees based on the variant allele frequency are proposed. Those methods can automatically produce the evolutionary trees. However, interpretation of the large number of trees from those methods remain unresolved. We propose the classification method for the cancer sub-clonal evolutionary trees. In this method, we adopt the tree space theory to analyse the tree set. The proposed method can identify the cluster structure based on the tree topology and edge length attributes. For the interpretation of the clusters, we also provide the method for evaluating the diversity of trees in the clusters, which gives an insight for the acceleration of sub-clonal expansion.

## Introduction

PhyC (Phylogenetic tree Clustering) is designed for classifying cancer evolutionary trees. 

The main inputs are
* list of edge matrix that must include the root
* list of edge length vector

and the main outpus are
* cluster assignments
* registered trees
* accerelation of sub-clonal expansions

In the input, the edge matrix represents cancer evolutionary topology and the edge length represents the number of additional Somatic Single Nucleotide Variants (SSNVs) from the parental sub-clone. Here is an illustration (Figure 1).

<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/ssnv.png" width="600" height="400" />

<strong>Figure 1.</strong> The model of cancer evolutionary trees

Our model assumes cancer clonal theory (Nowell,1976; Nik-Zainal et al., 2012) that
 
1. no mutation occurs twice in the course of cancer evolution
2. no mutation is ever lost.

Our classification method takes the inputs from those reconstruction methods below that outputs the various cancer evolutionary trees from the multiregional cancer sequencing per patiant and the excellent overview of those methods are described in (Beerenwinkel et al., 2014); <em>PhyloSub</em> (Jiao et al., 2013), <em>PyClone</em> (Roth et al., 2014), <em>SciClone</em> (Miller et al., 2014), <em>Clomial</em> (Zare et al., 2014), <em>Trap</em> (Strino et al., 2013), <em>SubcloneSeeker</em> (Qiao et al., 2014), and <em>LICHeE</em> (Popic et al., 2015) etc.

PhyC classifies those evolutinary trees based on the tree toplogies and edge length attributes under the framework of <em>tree space</em> (Birella,et al. 2001). Here is an example of a result (Figure 2).

<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/phyCMD.png" width="600" height="400" />

<strong>Figure 2.</strong> Example of the clustering

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

<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/phyC.plot.png" width="600" height="400" />

<strong>Figure 5.</strong> Example of non-hierarchical clustering

####diversity
The next function <em>diversity</em> is for interpretating clusters. We need some quantifications of trees in the clusters to evaluate them efficiently. We define the diversity as the number of sub-clones at a levels of SSNVs accumulations. In Figure 6, <em>h</em> indicates the normalized number of additional SSNVs (ranged [0,1]) and <em>g(h)</em> indicates the the number of subclones. This gives us the insights for accerelation of sub-clonal expansions. 

<img align="center" src="https://github.com/ymatts/PhyC/blob/master/img/diversity.png" width="600" height="400" />

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
