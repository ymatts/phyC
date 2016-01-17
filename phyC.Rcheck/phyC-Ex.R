pkgname <- "phyC"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "phyC-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('phyC')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("diversity")
### * diversity

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("diversity", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("phyC")
### * phyC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phyC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("phyCMD")
### * phyCMD

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phyCMD", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.phyC")
### * plot.phyC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.phyC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
