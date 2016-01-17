#' Internal function
#' 
#' @name get.index
#' @author Yusuke Matsui & Teppei Shimamura
#' @export

get.index <- function(edge){
  node <- unique(edge[,1])
  tips <- setdiff(as.vector(edge),node)
  root=rootFind(edge)
  node <- node[!node==root]
  res <- list(root = root,node=node,tips=tips)
  res
}
