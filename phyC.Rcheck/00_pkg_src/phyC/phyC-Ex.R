pkgname <- "phyC"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('phyC')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("diversity")
### * diversity

flush(stderr()); flush(stdout())

### Name: diversity
### Title: Evaluating diversity and plot
### Aliases: diversity

### ** Examples

library(phyC)
##generate edgeList and edgeLenList##
trees <- c(rmtree(5,3),rmtree(5,4))
edgeList <- lapply(trees,function(x)x$edge)
edgeLenList <- lapply(trees,function(x)x$edge.length)
##adopting phyC##
res <- phyC(edgeList,edgeLenList,cluster=2,type='nh')
div <- diversity(res)



cleanEx()
nameEx("phyC")
### * phyC

flush(stderr()); flush(stdout())

### Name: phyC
### Title: Clustering cancer evolutionary trees
### Aliases: phyC

### ** Examples

library(phyC)
##generate edgeList and edgeLenList##
trees <- c(rmtree(5,3),rmtree(5,4))
edgeList <- lapply(trees,function(x)x$edge)
edgeLenList <- lapply(trees,function(x)x$edge.length)
##adopting phyC##
res <- phyC(edgeList,edgeLenList,cluster=2,type='nh')



cleanEx()
nameEx("phyCMD")
### * phyCMD

flush(stderr()); flush(stdout())

### Name: phyCMD
### Title: Constructing configuration of trees in clusters and plot
### Aliases: phyCMD

### ** Examples

library(phyC)
##generate edgeList and edgeLenList##
trees <- c(rmtree(5,3),rmtree(5,4))
edgeList <- lapply(trees,function(x)x$edge)
edgeLenList <- lapply(trees,function(x)x$edge.length)
##adopting phyC##
res <- phyC(edgeList,edgeLenList,cluster=2,type='nh')
phyCMD(res)



cleanEx()
nameEx("plot.phyC")
### * plot.phyC

flush(stderr()); flush(stdout())

### Name: plot.phyC
### Title: Plot trees in clusters
### Aliases: plot.phyC

### ** Examples

library(phyC)
##generate edgeList and edgeLenList##
trees <- c(rmtree(5,3),rmtree(5,4))
edgeList <- lapply(trees,function(x)x$edge)
edgeLenList <- lapply(trees,function(x)x$edge.length)
##adopting phyC##
res <- phyC(edgeList,edgeLenList,cluster=2,type='nh')
plot.phyC(res)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
