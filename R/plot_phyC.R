phyC.plot <- function (obj, color = NULL, label = NULL, nPanel = c(5,5), ave = F,textsize=1,type = "unrooted",edge.width = 6) 
{
  if(ave){
    resolve_tree <- obj$ave_tree
  }else{
    resolve_tree <- obj$trees
  }
  cluster <- obj$cluster
  para <- par(list(mfrow = nPanel, mar = c(0, 0, 3, 0), bg = "#f0f0f0"))
  if (is.null(color)) {
    cn <- length(unique(cluster))
    if (cn < 3) {
      cn <- 3
    }
    color <- brewer.pal(cn, name = "Set2")
  }
  if (is.null(label)) {
    label <- paste0("Tree ", seq_along(resolve_tree))
  }
  for (i in seq_along(resolve_tree)) {
    plot.phylo(resolve_tree[[i]], type = type, edge.width = edge.width, 
               edge.color = color[cluster[i]], show.tip.label = F, 
               adj = 0.5)
    mtext(label[i], side = 3, line = 1, cex = textsize)
    idx <- get.index(resolve_tree[[i]]$edge)
    zerolen <- resolve_tree[[i]]$edge.length == 0
    zeronode <- resolve_tree[[i]]$edge[zerolen, 2]
    showtips <- idx$tips[is.na(match(idx$tips, zeronode))]
    tiplabels(tip = showtips, frame = "none", cex = 3, pch = 8, 
              col = "#66a61e")
    nodelabels(text = "", node = idx$root, frame = "none", 
               cex = 3, pch = 18, col = "black")
    box(which = "figure", lty = "solid", col = "white")
  }
  par(para)
}
