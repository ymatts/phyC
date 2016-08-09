#' Constructing phylogenetic trees.
#'
#' @name phyC
#' @param vaf Matrix of vaf profie. Each row and column must represent VAFs of genes and sampled region, respectively, and first column must normal cell.
#' @param thr Threshold of VAF, which is used for binarizing the VAF profile. If VAF>=thr then it is set to 1; otherwise 0.
#' @return tree Constructed phylogenetic trees.
#' @return edge Edge matrix of the constructed tree.
#' @return edge.length Edge length of the constructed tree.
#' @return node.label Labels of each noder representing sub-clones.
#' @return profile Mutation status of each sub-clone.
#' @details This function constructsparsimony trees such that nodes and edge length correspond to sub-clones and the number of SSNVs in the sub-clones. acctran in phangorn package is used to construct the parsimony trees.
#' @examples
#' library(phyC)
#' data(ccRCC)
#' trees <- par.tree(ccRCC)
#' res <- phyC(evol$edgeList,evol$edgeLenList,cluster=4,type='h',method="ward")
#' res$cluster
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#'
par.tree <- function(vaf,thr = 0.05){
  if("normal"%in%colnames(vaf)){
    ind <- which(colnames(vaf)=="normal")
    colnames(vaf)[ind] <- "Normal"
    }
  if(!all(apply(vaf,2,class)=="numeric")){
    stop("Invalid data type. Input data must be only vaf values\n")
  }
  vaf_bin <- vaf
  vaf_bin[vaf>=thr] <- 1
  vaf_bin[vaf<thr] <- 0
  vaf_bin <- as.data.frame(vaf_bin)
  phydat <- phyDat(vaf_bin,type="USER",levels=c(0,1))
  njtree <- NJ(dist.hamming(phydat))
  partree <- pratchet(phydat,trace = F)
  partree <- acctran(partree,phydat)
  tree <- as.phylo(partree)
  plot.phylo(tree,main=parname[i],type = "unrooted",direction = "rightwards",edge.width = 3,cex = 1.2)
  #tiplabels()
  
  
  g_undir <- as.igraph(tree,directed = F)
  tips <- tree$tip.label
  node_name <- names(V(g_undir))
  
  node <- seq_along(tips)
  normal <-  match(tips[match("Normal",tips)],names(V(g_undir)))
  subclone <- match(tips[-match("Normal",tips)],names(V(g_undir)))
  in_node <- (1:length(node_name))[-c(normal,subclone)]
  vpath <- sapply(subclone,function(x)get.shortest.paths(g_undir,normal,x)$vpath)
  
  edge_list <- vector("list",length(vpath))
  for(i in seq_along(vpath)){
    current_path <- as.numeric(vpath[[i]])
    edge <- matrix(NA,nrow=length(current_path)-1,ncol=2)
    for(j in 1:nrow(edge)){
      edge[j,] <- current_path[c(j,j+1)]
    }
    edge_list[[i]] <- edge
  }
  edge <- unique(do.call(rbind,edge_list))
  g_dir <- graph.edgelist(edge,directed=T)
  vertex.attributes(g_dir)$name <- rep("",length(node_name))
  vertex.attributes(g_dir)$name[c(normal,subclone)] <- node_name[c(normal,subclone)]
  vertex.attributes(g_dir)$name[-c(normal,subclone)] <- node_name[-c(normal,subclone)]
  #plot(g_dir, layout = layout.reingold.tilford(g_dir, root = normal))
  node_label <- vertex.attributes(g_dir)$name
  
  
  depth <- sapply(seq_along(node_name),function(x)length(get.shortest.paths(g_dir,normal,x)$vpath[[1]]))
  in_node_depth <- depth[in_node]
  unique_depth <- sort(unique(in_node_depth),decreasing=T)
  node_prof <- vector("list",length = length(node_name))
  names(node_prof) <- node_name
  r_id <- match(node_name[c(normal,subclone)],colnames(vaf_bin))
  node_prof[c(normal,subclone)] <- lapply(r_id,function(x)vaf_bin[,x])
  for(i in seq_along(unique_depth)){
    select_id <- which(unique_depth[i]==in_node_depth)
    for(j in seq_along(select_id)){
      child <- as.numeric(ego(g_dir,order=1,mode="out",nodes = in_node[select_id[j]])[[1]][-1])
      tmp_prof <- node_prof[[child[1]]]
      for(k in seq_along(child)[-1]){
        tmp_prof <- tmp_prof & node_prof[[child[k]]]
      }
      tmp_prof <- as.numeric(tmp_prof)
      ind <- in_node[select_id][j]
      node_prof[[ind]] <- tmp_prof
    }
  }
  
  edge_length <- rep(0,nrow(edge))
  for(i in 1:nrow(edge)){
    v1 <- node_name[edge[i,1]]
    v2 <- node_name[edge[i,2]]
    id1 <- match(v1,node_name)
    id2 <- match(v2,node_name)
    edge_length[i] <- sum(node_prof[[id2]]) - sum(node_prof[[id1]])
  }
  #edge_length <- replace(edge_length,edge_length!=0,1)
  out_tree <- list(tree=tree,edge=edge,edge.length=edge_length,node.label = node_label,profile=node_prof)
  out_tree #output
}
