#' Evaluating diversity and plot
#' 
#' @name diversity
#' @param obj Object resulted from phyC. 
#' @param color If plotit=TRUE, color parameter of each cluster in the plot.
#' @param plotit Whether plot the diversity or not (logical).
#' @return ind.div Diversity of each tree.
#' @return div Diversity of each cluster.
#' @examples
#' library(phyC)
#' data(ccRCC)
#' vaf <- lapply(vaf,function(x)x[,-(1:3)])
#' trees <- par.tree(vaf)
#' edgeList <- lapply(trees,function(x)x$edge)
#' edgeLenList <- lapply(trees,function(x)x$edge.length)
#' res <- phyC(edgeList,edgeLenList,cluster = 3,type = "h",method = "ward")
#' div <- diversity(res)
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#' 
diversity <- function(obj,color=NULL,label=NULL,plotit = T){
  edgeList <- obj$edgeList
  edgeLenList <- obj$edgeLenList
  cluster <- obj$cluster

if(is.null(color)){
      cn <- length(unique(cluster))
      if(cn<3){cn <- 3}
      color <- brewer.pal(cn,name = "Set2")
  }

  if(length(edgeList)!=length(edgeLenList)){print("the number of edgeList and edgLenList is different.");stop()}
  g.cumlen <- vector("list",length(length(edgeList)))
  for(i in seq_along(edgeList)){
    edge <- edgeList[[i]]
    edgelen <- edgeLenList[[i]]
    edgelen <- edgelen / sum(edgelen)
    g <- graph.edgelist(edge)
    g.idx <- get.index(edge)
    g.spath <- sapply(g.idx$tips,function(x)igraph::shortest_paths(g,from = g.idx$root,to = x,mode = "out")$vpath)
    cumlen <- matrix(0,nrow=length(unlist(g.spath)),ncol=2)
    cnt <- 1
    for(j in seq_along(g.spath)){
      v <- as.numeric(g.spath[[j]])
      ind <- get.edgeIndex(ref = edge,vec = g.spath[[j]])
      len <- c(0,cumsum(edgelen[ind]))
      for(k in seq_along(v)){
        cumlen[cnt,] <- c(v[k],len[k])
        cnt <- cnt + 1
      }
    }
    g.cumlen[[i]] <- unique(cumlen)
  }
  
  ax <- seq(0,1,0.01)
  nlist <- vector("list",length(g.cumlen))
  for(i in seq_along(g.cumlen)){
    cumlen <- g.cumlen[[i]]
    n <- rep(0,nrow(cumlen))
    for(j in seq_along(ax)){
      n[j] <- length(cumlen[cumlen[,2]<=ax[j],1])
    }
    nlist[[i]] <- n
  }
  
  div.mean <- vector("list",length(unique(cluster)))
  for(i in seq_along(div.mean)){
    temp <- nlist[cluster==i]
    temp <- do.call(rbind,temp)
    div.mean[[i]] <- apply(temp,2,mean)
  }
  
  if(plotit){
#     dat <- vector("list",length(nlist))
#     for(i in seq_along(nlist)){
#       tmp <- data.frame(x = ax,y = nlist[[i]],z = label[i],cluster = cluster[i])
#       dat[[i]] <- tmp
#     }
#     dat <- do.call(rbind,dat)
#     dat$cluster <- factor(dat$cluster,ordered = T)
# #    dat$z <- factor(dat$z)
#     
#     p <- ggplot(dat,aes(x,y,group = z,color=cluster)) + theme(legend.position="none") + scale_color_manual(values = color) + geom_line(size=1) + labs(list(x="The relative number of SSNVs accumulated",y="The number of sub-clones"))
#     #+coord_cartesian(xlim = c(min(dat$x), max(dat$x) + 0.1))
#     
#     #p <- p + geom_text_repel(data = subset(dat,x==max(x)),aes(label=z),size=3, nudge_x=45)
#     
#     
#     mean_dat <- vector("list",length(div.mean))
#     for(i in seq_along(div.mean)){
#       mean_dat[[i]] <- data.frame(x = ax, y = div.mean[[i]], cluster = i)
#     }
#     mean_dat <- do.call(rbind,mean_dat)
#     mean_dat$cluster <- factor(mean_dat$cluster,ordered = T)
#     p <- p + geom_line(aes(x,y,group = cluster),linetype= 2,data = mean_dat,size = 1.5)

    
    if(is.null(color)){color <- (min(cluster)):(max(cluster))}
    par(mar=c(4,4,1,1))
    maxnum <- max(unlist(nlist))+1
    plot(ax,nlist[[1]],type="n",xlim=c(0,1),ylim=c(0,maxnum),ann=F,axes=F)

    u <- par("usr") # The coordinates of the plot area
    rect(u[1], u[3], u[2], u[4], col="#f0f0f0", border=NA)
    
    grid(col = "white",lwd = 2,lty="solid")
    for(i in seq_along(nlist)){
      lines(ax,nlist[[i]],col=color[cluster[i]],lwd=1)
    }
    for(i in seq_along(div.mean)){
      lines(ax,div.mean[[i]],col=color[i],lwd=5,lty=2)
    }
    mtext("The relative number of SSNVs accumulated",side=1,line = 3,col="#969696",cex = 1.2)
    mtext("The number of sub-clones",side=2,line = 3,col = "#969696",cex = 1.2)
    axis(side = 1,col = "grey",font=5,col.ticks = "#969696")
    axis(side = 2,col = "grey",font=5,col.ticks = "#969696")
    
  }
  return(list(ind.div=nlist,div=div.mean))
}

