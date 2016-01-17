#' Internal function
#' 
#' @name create.meta
#' @author Yusuke Matsui & Teppei Shimamura
#' @export

create.metaTree <- function(depth){
  ntips <- 2^(depth-1) #the number of leaves
  tip.meta <- 1:ntips #index of leaf
  Nnode <- ntips - 1 #the number of internal nodes
  Nedge <- sum(2^(1:(depth-1))) #the number of edges
  g0 <- graph.tree(n = Nnode + ntips,children = 2,mode = "out")
  
  g.edge <- get.edgelist(g0)
  g.node_before <- unique(g.edge[,1])
  g.tips_before <- setdiff(unique(as.vector(g.edge)),g.node_before)
  g.root <- min(g.node_before)
  g <- graph.edgelist(g.edge)
  spath <- get.shortest.paths(g,from = g.root,to = g.tips_before,mode = "out")$vpath
  
  tips_after <- 1:ntips
  root_after <- ntips + 1
  node_after <- root_after:max(g.edge)
  
  temp.pair <- as.list(NULL)
  for(j in 1:length(spath)){
    diff.path <- sapply(spath,function(x)setdiff(union(spath[[j]],x),intersect(spath[[j]],x)))
    diff.len <- sapply(diff.path,length)
    pair.ind <- which(diff.len==2)
    if(length(pair.ind)==0){
      temp.pair <- c(temp.pair,list(as.vector(spath[[j]])))
    }else{
      temp.pair <- c(temp.pair,list(rbind(spath[[j]],spath[[pair.ind]])))
    }
  }
  
  tips_after <- 1:ntips
  root_after <- ntips + 1
  node_after <- root_after:max(g.edge)
  
  cnt.tips <- 1
  cnt.node <- 1
  cnt <- 1
  trans <- matrix(0,nrow = length(unique(as.vector(g.edge))), ncol = 2)
  pair <- temp.pair
  
  nodepath <- vector("list",length(pair))
  tips <- vector("list",length(pair))
  for(j in seq_along(pair)){
    if(class(pair[[j]])=="matrix"){
      temp.node <- pair[[j]][1,-ncol(pair[[j]])]
      temp.tips <- pair[[j]][,ncol(pair[[j]])]
    }else{
      temp.node <- pair[[j]][-length(pair[[j]])]
      temp.tips <- pair[[j]][length(pair[[j]])]
    }
    nodepath[[j]] <- temp.node
    tips[[j]] <- temp.tips
    
    for(k in seq_along(temp.node)){
      if(temp.node[k]%in%trans[,1]){next}
      trans[cnt,1] <- temp.node[k]
      trans[cnt,2] <- node_after[cnt.node]
      cnt.node <- cnt.node + 1
      cnt <- cnt + 1
    }
    for(k in seq_along(temp.tips)){
      if(temp.tips[k]%in%trans[,1]){next}
      trans[cnt,1] <- temp.tips[k]
      trans[cnt,2] <- tips_after[cnt.tips]
      cnt.tips <- cnt.tips + 1
      cnt <- cnt + 1
    }
  }
  
  new.g.edge <- g.edge
  for(j in 1:nrow(trans)){
    ind <- which(trans[j,1]==g.edge,arr.ind = T)
    new.g.edge[ind] <- trans[j,2]
  }
  res <- as.list(NULL)
  res$edge <- new.g.edge
  res$tip.label <- paste0("t",tip.meta)
  res$edge.length <- rep(0,times=nrow(res$edge))
  res$Nnode <- Nnode
  class(res) <- "phylo"
  res
}

#return index of root, internal nodes, leaves.
get.index <- function(edge){
  node <- unique(edge[,1])
  tips <- setdiff(as.vector(edge),node)
  root=min(node)
  node <- node[-which.min(node)]
  res <- list(root = root,node=node,tips=tips)
  res
}

#return index of edge length in phylo class
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
