#' Plot trees in clusters
#' 
#' @name plot.phyC
#' @param obj Object resulted from phyC. 
#' @param color Vector of color parameter of each cluster in the plot.
#' @param label Vector of labels of trees. Default is "Tree i"(i=1,2,...).
#' @examples
#' library(phyC)
#' ##generate edgeList and edgeLenList##
#' trees <- c(rmtree(5,3),rmtree(5,4))
#' edgeList <- lapply(trees,function(x)x$edge)
#' edgeLenList <- lapply(trees,function(x)x$edge.length)
#' ##adopting phyC##
#' res <- phyC(edgeList,edgeLenList,cluster=2,type='nh')
#' plot.phyC(res)
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#'
plot.phyC <- function(obj,color=NULL,label=NULL){
  
  if(length(resolve_tree)>=25){
    n <- length(resolve_tree)%%25
    cat("there are too many trees to plot them in one panel, then plots are divided into ",n," plots")
  }
  para <- par(list(mfrow=c(5,5),mar=c(1,1,3,1),bg="#f0f0f0"))
  
  resolve_tree <- obj$trees
  cluster <- obj$cluster

  if(is.null(color)){
    cn <- length(unique(cluster))
    if(cn<3){cn <- 3}
    color <- brewer.pal(cn,name = "Set2")
  }
  if(is.null(label)){label <- paste0("Tree ",seq_along(resolve_tree))}

  for(i in seq_along(resolve_tree)){  
    plot.phylo(resolve_tree[[i]],type = "unrooted",edge.width = 6,edge.color = color[cluster[i]],show.tip.label = F)
    mtext(label[i],side = 3,line = 1,cex=1)
    idx <- get.index(resolve_tree[[i]]$edge)    
    zerolen <- resolve_tree[[i]]$edge.length==0
    zeronode <- resolve_tree[[i]]$edge[zerolen,2]
    showtips <- idx$tips[is.na(match(idx$tips,zeronode))]            
    tiplabels(tip = showtips,frame = "none",cex=3,pch=8,col = "#66a61e")
    nodelabels(text = "",node = idx$root,frame = "none",cex=3,pch=18,col = "black")
    box(which="figure",lty="solid",col="white")
  }
  par(para)  
}
