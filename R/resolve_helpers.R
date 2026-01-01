#' Internal helper functions for resolving mono and multifurcations
#'
#' @name resolve_helpers
#' @author Yusuke Matsui & Teppei Shimamura
#' @keywords internal

#' Resolve mono-furcations by adding dummy nodes
#'
#' @param edge Edge matrix with edge lengths in third column
#' @return Updated edge matrix with mono-furcations resolved
#' @keywords internal
resolve_monofurcation <- function(edge) {
  tb <- table(edge[, 1])
  mono <- as.numeric(names(tb[tb == 1]))
  addnode <- -1

  for (j in seq_along(mono)) {
    edge <- rbind(edge, c(mono[j], addnode, 0))
    addnode <- addnode - 1
  }

  list(edge = edge, addnode = addnode)
}

#' Resolve multi-furcations by creating random binary subtrees
#'
#' @param edge Edge matrix with edge lengths in third column
#' @param addnode Current addnode counter
#' @param seed Random seed for reproducibility (default: 10)
#' @return Updated edge matrix with multi-furcations resolved
#' @keywords internal
resolve_multifurcation <- function(edge, addnode, seed = 10) {
  tb <- table(edge[, 1])
  multi <- as.numeric(names(tb[tb > 2]))

  for (j in seq_along(multi)) {
    edge_ind <- which(edge[, 1] == multi[j])
    tip.label <- edge[edge_ind, 2]
    len <- length(edge_ind)

    set.seed(seed)
    subtree <- ape::rtree(len, tip.label = tip.label)$edge
    subnode <- unique(subtree[, 1])
    subtree0 <- subtree

    # Renumber internal nodes
    for (k in 2:length(subnode)) {
      subtree[subtree0 == subnode[k]] <- addnode
      addnode <- addnode - 1
    }

    # Replace root of subtree with original multi-furcation node
    subst.ind <- subtree0 == subnode[1]
    subtree[subtree0 == subnode[1]] <- multi[j]
    temp.sub <- subtree0

    # Map tip labels back
    for (k in 1:length(tip.label)) {
      if (k == multi[j]) {
        subtree[temp.sub == k & !subst.ind] <- tip.label[k]
      } else {
        subtree[temp.sub == k] <- tip.label[k]
      }
    }

    # Add edge lengths
    subtree <- cbind(subtree, 0)
    temp <- edge[edge_ind, ]
    for (k in 1:length(tip.label)) {
      temp_len <- temp[temp[, 2] == tip.label[k], 3]
      subtree[subtree[, 2] == tip.label[k], 3] <- temp_len
    }

    edge <- edge[-edge_ind, ]
    edge <- rbind(edge, subtree)
  }

  edge
}

#' Build pair list from edge list for tree registration
#'
#' @param edge Edge matrix (without edge lengths)
#' @return List of path pairs
#' @keywords internal
build_pair_list <- function(edge) {
  idx <- get.index(edge)
  minnode <- min(edge)
  g.edge <- edge - minnode + 1
  g.idx <- get.index(g.edge)
  root.new <- g.idx$root
  g.tips_before <- g.idx$tips

  g <- igraph::graph.edgelist(g.edge)
  spath <- igraph::get.shortest.paths(g, from = root.new,
                                      to = g.tips_before, mode = "out")$vpath

  temp.pair <- as.list(NULL)
  for (j in seq_along(spath)) {
    diff.path <- sapply(spath, function(x)
      setdiff(union(spath[[j]], x), intersect(spath[[j]], x)))
    diff.len <- sapply(diff.path, length)
    pair.ind <- which(diff.len == 2)

    if (length(pair.ind) == 0) {
      temp.pair <- c(temp.pair, list(as.vector(spath[[j]])))
    } else {
      temp.pair <- c(temp.pair, list(rbind(spath[[j]], spath[[pair.ind]])))
    }
  }

  temp.pair
}

#' Build transformation matrix for node relabeling
#'
#' @param edge Edge matrix (without edge lengths)
#' @param pair Pair list from build_pair_list
#' @return Transformation matrix
#' @keywords internal
build_transformation <- function(edge, pair) {
  minnode <- min(edge)
  g.edge <- edge - minnode + 1
  g.idx <- get.index(g.edge)
  g.root <- g.idx$root
  g.tips_before <- g.idx$tips

  g <- igraph::graph.edgelist(g.edge)
  ntips <- length(g.tips_before)
  spath <- igraph::get.shortest.paths(g, from = g.root,
                                      to = g.tips_before, mode = "out")$vpath

  tips_after <- 1:ntips
  root_after <- ntips + 1
  node_after <- root_after:max(g.edge)

  cnt.tips <- 1
  cnt.node <- 1
  cnt <- 1
  trans <- matrix(0, nrow = length(unique(as.vector(g.edge))), ncol = 2)

  for (j in seq_along(pair)) {
    if (inherits(pair[[j]], "matrix")) {
      temp.node <- pair[[j]][1, -ncol(pair[[j]])]
      temp.tips <- pair[[j]][, ncol(pair[[j]])]
    } else {
      temp.node <- pair[[j]][-length(pair[[j]])]
      temp.tips <- pair[[j]][length(pair[[j]])]
    }

    for (k in seq_along(temp.node)) {
      if (temp.node[k] %in% trans[, 1]) next
      trans[cnt, 1] <- temp.node[k]
      trans[cnt, 2] <- node_after[cnt.node]
      cnt.node <- cnt.node + 1
      cnt <- cnt + 1
    }

    for (k in seq_along(temp.tips)) {
      if (temp.tips[k] %in% trans[, 1]) next
      trans[cnt, 1] <- temp.tips[k]
      trans[cnt, 2] <- tips_after[cnt.tips]
      cnt.tips <- cnt.tips + 1
      cnt <- cnt + 1
    }
  }

  trans
}

#' Convert edge matrix to phylo object
#'
#' @param edge Edge matrix with edge lengths in third column
#' @return phylo object
#' @keywords internal
edge_to_phylo <- function(edge) {
  edge_only <- edge[, -3, drop = FALSE]
  node <- unique(edge_only[, 1])
  tip.label <- setdiff(as.vector(edge_only), node)
  edge.length <- edge[, 3]
  Nnode <- length(node)

  tree <- list(
    edge = edge_only,
    tip.label = tip.label,
    edge.length = edge.length,
    Nnode = Nnode
  )
  class(tree) <- "phylo"

  tree
}
