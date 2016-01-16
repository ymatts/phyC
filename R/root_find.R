#' Internal function
#' 
#' @name resolve.monomulti
#' @docType package
#' @import igraph,ape
#' @author Yusuke Matsui & Teppei Shimamura
#' @export


rootFind <- function(edge){
  node <- unique(edge[,1])
  ind <- !node%in%edge[,2]
  if(sum(ind)==1){
    root <- node[ind]
  }else{
    print("there is not root node")
    stop()
  }
  root
}
