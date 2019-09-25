library(clusterProfiler)
library(paxtoolsr)
library(org.Hs.eg.db)
library(tidygraph)
library(tidytext)

library(readxl)

source("modified_clusterProfiler_functions.R")

# Return variable as string 
# From: https://stackoverflow.com/questions/14577412/how-to-convert-variable-object-name-into-string
#curVar <- function(x) { deparse(substitute(x)) }

# Load Data 
data(stop_words)
data(geneList, package="DOSE")
dose_geneList <- geneList

up_erbb2_lo <- read_xlsx("breast_data/nature19328-s5.xlsx", skip = 2, sheet = "Up in ERBB2 low")
up_erbb2_hi <- read_xlsx("breast_data/nature19328-s5.xlsx", skip = 2, sheet = "Up in ERBB2 high")

geneList <- up_erbb2_lo$Symbol

# geneList <- sample(geneList, 5000)

# Parameters 
maxTags <- 5

# Convert Pathway Commons GMT to ClusterProfiler format 
tmpFile <- "PathwayCommons11.All.hgnc.gmt.gz"
if(!file.exists(tmpFile)) {
  download.file("https://www.pathwaycommons.org/archives/PC2/v11/PathwayCommons11.All.hgnc.gmt.gz", tmpFile)
} 
gmt <- readGmt(tmpFile, returnInfo = TRUE)  
#gmt <- gmt[1:2000]

geneSetList <- lapply(seq_along(gmt), function(x, n, i) { 
  tmp <- x[[i]]
  data.frame(id=n[i], name=tmp[["name"]], gene=tmp[["geneSet"]], stringsAsFactors=FALSE)
}, x=gmt, n=names(gmt))

tmp <- do.call("rbind", geneSetList)
rownames(tmp) <- 1:nrow(tmp) # For convenience 

pc2gene <- tmp[, c("id", "gene")]
pc2name <- tmp[, c("id", "name")]

# ENRICHMENT MAP ----
## For DOSE 
#tmpDeg <- names(geneList)[abs(geneList) > 2]
#deg <- mapIds(org.Hs.eg.db, keys=tmpDeg, column="SYMBOL", keytype="ENTREZID", multiVals="first")
deg <- geneList
enrichOutput <- enricher(deg, pvalueCutoff=0.05, minGSSize=10, maxGSSize=500, TERM2GENE=pc2gene, TERM2NAME=pc2name)
a <- enrichOutput@result

em <- emapplot(enrichOutput, showCategory=20)
em

g <- emapplotModGraph(enrichOutput, showCategory=20)
#plot(g)

# Modify graph 
tmp <- components(g) 
gComponents <- data.frame(name=names(tmp$membership), component=as.numeric(tmp$membership))

#tmp <- cluster_fast_greedy(g) 
#gComponents <- data.frame(name=tmp$names, component=as.numeric(tmp$membership))

g1 <- as_tbl_graph(g)
nodes <- g1 %>% activate(nodes) %>% as_tibble
edges <- g1 %>% activate(edges) %>% as_tibble
nodeOrder <- nodes$name
limits <- c(min(nodes$color), max(nodes$color))

tmpNodes <- merge(nodes, gComponents, by="name", all.x=TRUE)
nodesOrdered <- tmpNodes[match(nodeOrder, tmpNodes$name),]

table(nodesOrdered$component)

g2 <- tbl_graph(
  nodes = nodesOrdered, 
  edges = edges, 
  directed = FALSE
)

tags <- sapply(1:max(nodesOrdered$component), function(i) {
  #i <- 2
  tmp <- nodesOrdered[nodesOrdered$component == i, "name"]
  
  tmp <- gsub("[[:punct:]]", " ", tmp)
  tmp <- gsub("  ", " ", tmp)
  
  a <- strsplit(tmp, " ") %>% unlist %>% unique 
  b <- c(tmp, a) %>% unique
  b <- a %>% unique
  b <- b[!(b %in% stop_words$word)]
  b <- b[nchar(b) > 3]
  
  c <- sapply(b, function(x) {
    grepl(x, tmp) %>% which %>% length 
  }) %>% sort(decreasing = TRUE)
  
  if(max(c) > 1) {
    d <- c[c > 1] %>% sort(decreasing = TRUE)
  } else {
    d <- c %>% sort(decreasing = TRUE)
  }
  
  maxTags <- ifelse(length(d) > maxTags, maxTags, length(d))
  tags <- paste(names(d)[1:maxTags], collapse = ", ")
  tags
})

tags <- paste0("Tags: ", tags)

emapplotModPlot(g2, componentHighlight=1:3, tags=tags, limits=limits)
