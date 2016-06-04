#' Internal function
#' 
#' @name regis
#' @author Yusuke Matsui & Teppei Shimamura
#' @export

regis <- function(edgelist,edgeLenList){
  resolv <- resolve.monomulti(edgelist,edgeLenList)
  maxdep <- depth.max(resolv)
  cat("creating maximum tree\n")
  cat(paste0("depth=",maxdep," and #leaves=",2^(maxdep),"\n"))
  #meta <- create.metaTree(maxdep+1)
  meta <- stree(ntips <- 2^(maxdepth-1+1),"balanced")
  cat("starting encoding trees to maximum tree\n")
  regis <- vector("list",length(resolv))
  for(i in seq_along(resolv)){
    regis[[i]] <- meta.regis.tree(meta,resolv[[i]])
    cat(i,"th registration finished\n")
  }
  regis
}
