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

emapplotModPlot <- function(g, limits, color="p.adjust", layout="nicely", tags=NULL) {
  p <- ggraph(g, layout=layout)
  
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha=0.8, aes_(width=~I(width)), colour='darkgrey')
  }
  
  p <- p + geom_node_point(aes(fill=color, size=size), color="black", shape=21, stroke=1)
  #p <- p + geom_node_point(aes(fill=color, size=size), stroke=2)
  p <- p + theme_graph()
  p <- p + scale_fill_gradient2(low="cyan", mid="white", high="magenta", midpoint=0,
                         limits=limits,
                         breaks= c(min(limits), max(limits)),
                         labels= c("ERBB2 Lo", "ERBB2 Hi"),
                         name="Enrichment Score")
  p <- p + scale_size(range=c(3, 12), name="Pathway Size")
   
  #   #geom_node_text(aes_(label=~name), repel=TRUE) + 
  #   #geom_node_text(aes_(label=~component), repel=TRUE) + 
 
  if(!is.null(tags)) {
    p <- p + lapply(1:nrow(tags), function(i) {
      #cat("X: ", a, " Y: ", b, "\n")
      geom_mark_hull(aes(x, y, label=tags[i, "tag"], filter=component == tags[i, "component"]), 
                     con.cap = 0, concavity=20)
    })
  }
   
  return(p)
}