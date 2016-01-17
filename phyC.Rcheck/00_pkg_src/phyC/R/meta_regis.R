#' Internal function
#' 
#' @name meta.regis.tree
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#'  
meta.regis.tree <- function(meta,target){

  edge.meta <- meta$edge
  edge.target <- target$edge
  
  g.meta <- igraph::graph.edgelist(edge.meta)
  g.target <- igraph::graph.edgelist(edge.target)
  index.meta <- get.index(edge.meta)
  index.target <- get.index(edge.target) 
  spath.meta <- sapply(index.meta$tips,function(x)shortest_paths(g.meta,from = index.meta$root,to = x)$vpath)
  spath.target <- sapply(index.target$tips,function(x)shortest_paths(g.target,from = index.target$root,to = x)$vpath)
  maxdepth <- length(spath.meta[[1]])
  target.depth <- sapply(spath.target,length)
  target.length <- sapply(spath.target,function(x)sum(target$edge.length[get.edgeIndex(target$edge,x)]))
  target.profile <- cbind(target.depth,target.length)
  ord <- order(target.profile[,1],target.profile[,2],decreasing=T)
  trans <- matrix(0,nrow=max(as.vector(target$edge)),ncol=2)
  
  cnt <- 1
  tree.regis <- meta
  for(i in seq_along(spath.target)){
    temp.target <- spath.target[[ord[i]]]
    depth.target <- target.profile[ord[i],1]
    dgap <- maxdepth - depth.target
    if(dgap==0){print("not enough depth of meta-tree");stop()}
    if(i!=1){
      cnt <- cnt + 2^(maxdepth - depth.target)
    }
    temp.meta <- spath.meta[[cnt]][1:depth.target]
    idx.meta <- get.edgeIndex(meta$edge,temp.meta)
    idx.target <- get.edgeIndex(target$edge,temp.target)
    tree.regis$edge.length[idx.meta] <- target$edge.length[idx.target]
  }
  tree.regis
}
