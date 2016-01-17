#' Internal function
#' 
#' @name depth.max
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#' 
depth.max <- function(trees){
  edgelist <- lapply(trees,function(x)x$edge)
  g.idx <- lapply(edgelist,function(x)get.index(x))
  g <- lapply(edgelist,graph.edgelist)
  maxdepth <- rep(NA,length(trees))
  for(i in seq_along(trees)){
    temp <- get.shortest.paths(g[[i]],from = g.idx[[i]]$root,to = g.idx[[i]]$tips)$vpath
    maxdepth[i] <- max(sapply(temp,length))
  }
  max(maxdepth)
}
