#' Plot trees in clusters
#' 
#' @name phyC.plot
#' @param obj Object resulted from phyC. 
#' @param color Vector of color parameter of each cluster in the plot.
#' @param label Vector of labels of trees. Default is "Tree i"(i=1,2,...).
#' @param type Type of tree to be plotted.  It must be one of "unrooted"(the default), "phylogram", "cladogram", "fan", "radial".
#' @examples
#' library(PhyC)
#' data(evol)
#' res <- phyC(evol$edgeList,evol$edgeLenList,cluster=4,type='nh')
#' phyC.plot(res)
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#'
phyC.plot <- function(obj,color=NULL,label=NULL,type="unrooted"){
  resolve_tree <- obj$trees
#   normalize <- obj$normalize
#   if(normalize=="total"){
#     for(i in seq_along(resolve_tree)){
#       resolve_tree[[i]]$edge.length <- resolve_tree[[i]]$edge.length / sum(resolve_tree[[i]]$edge.length)
#     }
#   }else if(normalize=="subclone"){
#     for(i in seq_along(resolve_tree)){
#       resolve_tree[[i]]$edge.length <- resolve_tree[[i]]$edge.length / resolve_tree[[i]]$edge.length[1]
#     }
#   }else if(normalize=="trunc"){
#     for(i in seq_along(resolve_tree)){
#       resolve_tree[[i]]$edge.length <- resolve_tree[[i]]$edge.length / sum(resolve_tree[[i]]$edge.length[-1])
#     }
#   }
  
  cluster <- obj$cluster

  if(length(resolve_tree)>=25){
    n <- length(resolve_tree)%%25
    cat("there are too many trees to plot them in one panel, then plots are divided into ",n," plots")
  }
  #para <- par(list(mfrow=c(5,5),mar=c(0,0,3,0),bg="#f0f0f0"))
  para <- par(list(mfrow=c(5,5),mar=c(0,0,3,0),bg="#f0f0f0"))

  if(is.null(color)){
    cn <- length(unique(cluster))
    if(cn<3){cn <- 3}
    color <- brewer.pal(cn,name = "Set2")
  }
  if(is.null(label)){label <- paste0("Tree ",seq_along(resolve_tree))}

  for(i in seq_along(resolve_tree)){  
    plot.phylo(resolve_tree[[i]],type = type,edge.width = 6,edge.color = color[cluster[i]],show.tip.label = F,adj = 0.5)
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
