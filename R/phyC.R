#' Clustering cancer evolutionary trees
#' 
#' @name phyC
#' @docType package
#' @import igraph,ape
#' @param edgeList List of edge matrix that include the root.
#' @param edgeLenList List of edge length vector. Index of vector should be corresponding to row index of edge matrix.
#' @param cluster The number of clusters.
#' @param type Clustering type. Selecht 'nh' (non hierarchical) or 'h'(hierarchical with ward's method). 
#' @param edgeList List of edge matrix of input. 
#' @param edgeLenList List of edge length vector of input. 
#' @return trees Input trees.
#' @return regis.tree Registered trees.
#' @return cluster Index of the clusters.
#' @return dist Distance betwee trees.
#' @details This function perform the registration and the clustering. In the registration, we resolve the mono- and multi-furcation trees and we complete the number of leaves among the trees. In the resigration, identical tree toplogies are regarded identical even if the label are different. In the clustering, we classify the trees into predifined number of subsets. When choose type='nh', we perform the Ward's clustering.
#' @examples
#' library(phyC)
#' ##generate edgeList and edgeLenList##
#' trees <- c(rmtree(5,3),rmtree(5,4))
#' edgeList <- lapply(trees,function(x)x$edge)
#' edgeLenList <- lapply(trees,function(x)x$edge.length)
#' ##adopting phyC##
#' res <- phyC(edgeList,edgeLenList,cluster=2,type='nh')
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#' 
phyC <- function(edgeList,edgeLenlist,cluster,type='nh',method=NULL){
  library(igraph)
  library(ape)
  cat("staring to resolve mono and multifurcations \n")
  new.edgeList <- vector("list",length(edgeList))
  for(i in 1:length(edgeList)){
    edge <- cbind(edgeList[[i]],edgeLenList[[i]])
    idx <- get.index(edge[,-3,drop=F])
    root_before <- idx$root
    new.edge <- edge
    maxtip <- max(idx$tips)
    tips_before <- idx$tips
    tb <- table(edge[,1])
    mono <- as.numeric(names(tb[tb==1]))
    addnode <- -1
    for(j in seq_along(mono)){
      new.edge <- rbind(new.edge,c(mono[j],addnode,0))
      addnode <- addnode - 1
    }
    
    multi <- as.numeric(names(tb[tb > 2]))
    
    for(j in seq_along(multi)){
      edge_ind <- which(new.edge[,1]==multi[j])
      tip.label <- new.edge[edge_ind,2]
      len <- length(edge_ind)
      set.seed(10)
      subtree <- ape::rtree(len,tip.label=tip.label)$edge
      subnode <- unique(subtree[,1])
      
      subtree0 <- subtree
      for(k in 2:length(subnode)){
        subtree[subtree0==subnode[k]] <- addnode
        addnode <- addnode - 1
      }
      
      subst.ind <- subtree0==subnode[1]
      subtree[subtree0==subnode[1]] <- multi[j]
      temp.sub <- subtree0
      
      for(k in 1:length(tip.label)){
        if(k==multi[j]){
          subtree[temp.sub==k & !subst.ind] <- tip.label[k]
        }else{
          subtree[temp.sub==k] <- tip.label[k]
        }
      }
      
      subtree <- cbind(subtree,0)
      
      temp <- new.edge[edge_ind,]
      for(k in 1:length(tip.label)){
        temp_len <- temp[temp[,2] == tip.label[k],3]
        subtree[subtree[,2] == tip.label[k],3] <- temp_len
      }
      new.edge <- new.edge[-edge_ind,]
      new.edge <- rbind(new.edge,subtree)
    }
    new.edgeList[[i]] <- new.edge
  }
  
  pair.list <- vector("list",length(new.edgeList))
  for(i in seq_along(new.edgeList)){
    edge <- new.edgeList[[i]]
    idx <- get.index(edge[,-3,drop=F])
    tips_before <- idx$tips
    minnode <- min(edge[,1:2,drop=F])
    node_before <- idx$node
    g.edge <- edge[,-3,drop=F] - minnode + 1
    g.idx <- get.index(g.edge)
    root.new <- g.idx$root
    g.tips_before <- g.idx$tips
    g <- igraph::graph.edgelist(g.edge)
    ntips <- length(g.tips_before)
    spath <- igraph::get.shortest.paths(g,from = root.new,to = g.tips_before,mode = "out")$vpath
    #spath <- sapply(g.tips_before,function(x)shortest_paths(g,from = root.new,to = x,mode = "out")$vpath)
    
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
    pair.list[[i]] <- temp.pair
  }
  
  
  trans.list <- vector("list",length(pair.list))
  for(i in seq_along(pair.list)){
    edge <- new.edgeList[[i]]
    idx <- get.index(edge[,-3,drop=F])
    tips_before <- idx$tips
    root_before <- idx$root
    
    maxnode <- max(edge[,1:2,drop=F])
    minnode <- min(edge[,1:2,drop=F])
    
    g.edge <- edge[,-3,drop=F] - minnode + 1
    g.idx <- get.index(g.edge)
    g.root <- g.idx$root
    
    g.tips_before <- g.idx$tips
    g <- igraph::graph.edgelist(g.edge)
    ntips <- length(g.tips_before)
    spath <- igraph::get.shortest.paths(g,from = g.root,to = g.tips_before,mode = "out")$vpath
    
    tips_after <- 1:ntips
    root_after <- ntips + 1
    node_after <- root_after:max(g.edge)
    
    cnt.tips <- 1
    cnt.node <- 1
    cnt <- 1
    trans <- matrix(0,nrow = length(unique(as.vector(g.edge))), ncol = 2)
    pair <- pair.list[[i]]
    
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
    trans.list[[i]] <- trans
  }
  
  new.g.edge.list <- vector("list",length(new.edgeList))
  for(i in seq_along(new.edgeList)){
    edge <- new.edgeList[[i]]
    
    idx <- get.index(edge[,-3,drop=F])
    tips_before <- idx$tips
    root_before <- idx$root
    
    maxnode <- max(edge[,1:2,drop=F])
    minnode <- min(edge[,1:2,drop=F])
    
    g.edge <- edge[,-3,drop=F] - minnode + 1
    g.idx <- get.index(g.edge)
    g.root <- g.idx$root
    new.g.edge <- g.edge
    trans <- trans.list[[i]]
    for(j in 1:nrow(trans)){
      ind <- which(trans[j,1]==g.edge,arr.ind = T)
      new.g.edge[ind] <- trans[j,2]
    }
    new.g.edge.list[[i]] <- cbind(new.g.edge,edge[,3])
  }
  
  treelist <- vector("list",length(new.g.edge.list))
  for(i in seq_along(new.g.edge.list)){
    edge <- new.g.edge.list[[i]][,-3,drop=F]
    node <- unique(edge[,1])
    tip.label <- setdiff(as.vector(edge),node)
    edge.length <- new.g.edge.list[[i]][,3]
    Nnode <- length(node)
    temp.tree <- list(edge = edge,tip.label=tip.label,edge.length=edge.length,Nnode=Nnode)
    class(temp.tree) <- "phylo"
    treelist[[i]] <- temp.tree
  }
  
  for(i in seq_along(treelist)){
    treelist[[i]]$edge.length <- treelist[[i]]$edge.length / sum(treelist[[i]]$edge.length)
  }
  
  maxdep <- depth.max(treelist)
  cat("creating maximum tree\n")
  cat(paste0("depth=",maxdep," and #leaves=",2^(maxdep),"\n"))
  meta <- create.metaTree(maxdep+1)
  cat("starting encoding trees to maximum tree\n")
  regis <- vector("list",length(treelist))
  for(i in seq_along(treelist)){
    regis[[i]] <- meta.regis.tree(meta,treelist[[i]])
    cat(i,"th registration finished\n")
  }  
  coord <- t(sapply(regis,function(x)x$edge.length))
  d <- dist(coord)
  if(type=="nh"){
  cat("non-hierarchical clustering")  
  km <- kmeans(coord,cluster)
  cls <- km$cluster
 }else if(type=="h"){
  cat("hierarchical clustering")  
  if(method==NULL){method <- "ward.D2"}
  hc <- hclust(dist(coord),method=method)
  if(dendro==TRUE){plot(hc)}
  cls <- cutree(hc,clster)
 }
 return (list(trees=treelist,regis.trees=regis,cluster=cls,dist=d,edgeList=edgeList, edgeLenList=edgeLenList))
}
