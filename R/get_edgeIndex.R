#' Internal function
#' 
#' @name get.edgeIndex
#' @author Yusuke Matsui & Teppei Shimamura
#' @export

get.edgeIndex <- function(ref,vec){
  edge <- matrix(0,nrow=length(vec)-1,ncol=2)
  for(i in 1:(length(vec)-1)){
    edge[i,] <- c(vec[i],vec[i+1])
  }
  ind <- rep(0,length(nrow(edge)))
  for(i in 1:nrow(edge)){
    ind[i] <- which(apply(ref,1,function(x)all(x==edge[i,])))
  }
  ind
}
