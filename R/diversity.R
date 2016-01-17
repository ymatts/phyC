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
#' data(evol)
#' res <- phyC(edgeList,edgeLenList,cluster=2,type='nh')
#' div <- diversity(res)
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#' 
diversity <- function(obj,color=NULL,plotit = T){
  edgeList <- obj$edgeList
  edgeLenList <- obj$edgeLenList
  cluster <- obj$cluster
  
  if(length(edgeList)!=length(edgeLenList)){print("the number of edgeList and edgLenList is different.");stop()}
  g.cumlen <- vector("list",length(length(edgeList)))
  for(i in seq_along(edgeList)){
    edge <- edgeList[[i]]
    edgelen <- edgeLenList[[i]]
    edgelen <- edgelen / sum(edgelen)
    g <- graph.edgelist(edge)
    g.idx <- get.index(edge)
    g.spath <
      duver- sapply(g.idx$tips,function(x)igraph::shortest_paths(g,from = g.idx$root,to = x,mode = "out")$vpath)
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
    if(is.null(color)){color <- (min(cluster)):(max(cluster))}
    par(mar=c(4,4,1,1))
    maxnum <- max(unlist(nlist))+1
    plot(ax,nlist[[1]],type="n",xlim=c(0,1),ylim=c(1,maxnum),ann=F)
    grid()
    for(i in seq_along(nlist)){
      lines(ax,nlist[[i]],col=color[cluster[i]],lwd=1)
    }
    for(i in seq_along(div.mean)){
      lines(ax,div.mean[[i]],col=color[i],lwd=5,lty=2)
    }
    mtext("h",side=1,line = 3)
    mtext("g(h)",side=2,line = 3)
  }
  return(list(ind.div=nlist,div=div.mean))  
}

