#' Constructing configuration of trees in clusters and plot
#'
#' @name phyCMD
#' @param obj Object resulted from phyC.
#' @param color Vector of color parameter of each cluster in the plot.
#' @param label Vector of labels of trees. Default is "Tree i"(i=1,2,...).
#' @param img.width(/img.height) Image width (/ height) of the trees to be overlayed on the Euclidean space. The unit is "px". Default is 200.
#' @param size Size parameter of trees to be plotted in the Euclidean space. Usually, 0 < size < 1. Default is 0.15.
#' @param type Type of tree to be plotted.  It must be one of "unrooted"(the default), "phylogram", "cladogram", "fan", "radial".
#' @return dist Distance used for configuration
#' @return coord Coordinate of trees in the Euclidean space.
#' @details This function performs classical multidimensional scaling with tree distance. The resulting plot includes the trees overlayed on the Euclidean coordinates.
#' @examples
#' library(PhyC)
#' data(evol)
#' res <- phyC(evol$edgeList,evol$edgeLenList,cluster=4,type='nh')
#' phyCMD(res)
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#'

phyCMD <- function(obj,color=NULL,label=NULL,img.width=200,img.height=200,size=0.15,type="unrooted"){
      
      resolve_tree <- obj$trees
#       normalize <- obj$normalize
#       if(normalize=="total"){
#         for(i in seq_along(resolve_tree)){
#           resolve_tree[[i]]$edge.length <- resolve_tree[[i]]$edge.length / sum(resolve_tree[[i]]$edge.length)
#         }
#       }else if(normalize=="subclone"){
#         for(i in seq_along(resolve_tree)){
#           resolve_tree[[i]]$edge.length <- resolve_tree[[i]]$edge.length / resolve_tree[[i]]$edge.length[1]
#         }
#       }else if(normalize=="trunc"){
#         for(i in seq_along(resolve_tree)){
#           resolve_tree[[i]]$edge.length <- resolve_tree[[i]]$edge.length / sum(resolve_tree[[i]]$edge.length[-1])
#         }
#       }
      
      cluster <- obj$cluster
      d <- obj$dist
      
      if(is.null(color)){
          cn <- length(unique(cluster))
          if(cn<3){cn <- 3}
          color <- brewer.pal(cn,name = "Set2")
      }
      bgcol <- "#f0f0f0"
      
    if(is.null(label)){label <- paste0("Tree ",seq_along(resolve_tree))}
    
    outdir <- tempdir()
    output <- paste0(outdir,"phylo",seq_along(resolve_tree),".png")
    for(i in seq_along(resolve_tree)){
        png(filename = output[i],width = img.width,height = img.height,bg="transparent")
        para <- par(mar=c(0,0,0,0))
        plot.phylo(resolve_tree[[i]],type = type,edge.width = 6,edge.color = color[cluster[i]],show.tip.label = F,adj = 0.5)
        #mtext(label[i],side = 1,line = -3,cex=1)
        idx <- get.index(resolve_tree[[i]]$edge)
        zerolen <- resolve_tree[[i]]$edge.length==0
        zeronode <- resolve_tree[[i]]$edge[zerolen,2]
        showtips <- idx$tips[is.na(match(idx$tips,zeronode))]
        tiplabels(tip = showtips,frame = "none",cex=3,pch=8,col = "#66a61e")
        nodelabels(text = "",node = idx$root,frame = "none",cex=3,pch=18,col = "black")
        
        dev.off()
        par(para)
    }
    
    cmd <- cmdscale(d,k=2,eig=F)
    df <- as.data.frame(cmd)
    colnames(df) <- c("x1","x2")
    
    ranx <- range(as.vector(cmd[,1]))
    rany <- range(as.vector(cmd[,2]))
    
    p <- ggplot(df,aes(x1,x2,label=label),cmd) + geom_blank() +
    theme(panel.background = element_rect(fill = bgcol, color = bgcol, size = 2))+
    xlim(c(ranx[1]*1.2,ranx[2]*1.2))+
    ylim(c(rany[1]*1.2,rany[2]*1.2))
    
    
    radi <- size*min(abs(max(ranx[2],ranx[1])-min(ranx[2],ranx[1])),abs(max(rany[2],rany[1])-min(rany[2],rany[1])))
    for(i in seq_along(output)){
        img <- readPNG(output[i])
        p <- p + annotation_raster(img,xmin = cmd[i,1]-radi,xmax= cmd[i,1]+radi,ymin=cmd[i,2]-radi,ymax = cmd[i,2]+radi,interpolate = T)
    }
    p <- p + geom_text()
    plot(p)
    try(invisible(file.remove(output)))
    return(list(coord=cmd,dist=d))
}
