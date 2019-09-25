library(igraph)
library(ggraph)
library(DOSE)
library(reshape2)
library(ggforce)

# From clusterProfiler
update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  
  return(n)
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

emapplotModGraph <- function(x, showCategory = 30, color="p.adjust", layout = "kk", ...) {
  n <- update_n(x, showCategory)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  if (is.numeric(n)) {
    y <- y[1:n,]
  } else {
    y <- y[match(n, y$Description),]
    n <- length(n)
  }
  
  
  if (n == 0) {
    stop("no enriched term found...")
  } else if (n == 1) {
    g <- graph.empty(0, directed=FALSE)
    g <- add_vertices(g, nv = 1)
    V(g)$name <- y$Description
    V(g)$color <- "red"
    return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
  } else {
    id <- y[,1]
    geneSets <- geneSets[id]
    
    n <- nrow(y) #
    w <- matrix(NA, nrow=n, ncol=n)
    colnames(w) <- rownames(w) <- y$Description
    
    for (i in 1:n) {
      for (j in i:n) {
        w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }
    
    wd <- melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    wd <- wd[!is.na(wd[,3]),]
    g <- graph.data.frame(wd[,-3], directed=FALSE)
    E(g)$width=sqrt(wd[,3] * 5)
    g <- delete.edges(g, E(g)[wd[,3] < 0.2])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
    
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    
    colVar <- y[idx, color]
    V(g)$color <- colVar
  }

  return(g)
}

emapplotModPlot <- function(g, limits, color="p.adjust", layout="kk", componentHighlight=1, tags="tags") {
  p <- ggraph(g, layout=layout)
  
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
  }
  
  p + geom_node_point(aes_(color=~color, size=~size, stroke=1)) +
    
    #geom_node_text(aes_(label=~name), repel=TRUE) + 
    #geom_node_text(aes_(label=~component), repel=TRUE) + 
    geom_mark_hull(con.cap = 0, concavity=20, aes(x, y, label=tags[1], filter=component %in% componentHighlight[1])) +
    geom_mark_hull(con.cap = 0, concavity=20, aes(x, y, label=tags[2], filter=component %in% componentHighlight[2])) +
    geom_mark_hull(con.cap = 0, concavity=20, aes(x, y, label=tags[3], filter=component %in% componentHighlight[3])) +
    # geom_mark_hull(concavity=20, color="white", aes(x, y, label=tags) +
    theme_void() +
    scale_color_continuous(low="red", high="white", limits=limits, name=color, guide=guide_colorbar()) +
    ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range=c(3, 8))
}