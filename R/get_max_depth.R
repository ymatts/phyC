#' Internal function
#' 
#' @name get.maxdepth
#' @docType package
#' @import igraph,ape
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#' 
get.maxdepth <- function(trees){
  edgelist <- lapply(trees,function(x)x$edge)
  glist <- lapply(edgelist,function(x)igraph::graph.edgelist(x))
  idx.list<- lapply(edgelist,get.index)
  maxdep.list <- rep(NA,length(idx.list))
  for(i in seq_along(idx.list)){
    temp <- sapply(idx.list[[i]]$tips,function(x)igraph::shortest_paths(glist[[i]],from = idx.list[[i]]$root,x)$vpath)
    maxdep.list[i] <- max(sapply(temp,length))
  }
  maxdep.list
}
