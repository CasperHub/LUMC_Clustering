### Import required packages  
suppressPackageStartupMessages({
  if (!require("rlang")) install.packages("rlang")
  install.packages("dynamicTreeCut")
  library("dynamicTreeCut")
  if (!require("limma")) install.packages("limma")
  if (!require("edgeR")) install.packages("edgeR")
  if (!require("bluster")) install.packages("bluster")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("scRNAseq")
  BiocManager::install("scuttle")
  BiocManager::install("scater")
  BiocManager::install("scran")
  devtools::install_github("LUMC/dgeAnalysis")
  library("scRNAseq")
  library("scuttle")
  library("scater")
  library("scran")
  library("dgeAnalysis")
})

# Load data
sce <- ZeiselBrainData()

# Trusting the authors' quality control, and going straight to normalization.

sce <- logNormCounts(sce)

# Feature selection based on highly variable genes.

dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n=1000)

# Dimensionality reduction for work (PCA) and pleasure (t-SNE).
set.seed(1000)

sce <- runPCA(sce, ncomponents=20, subset_row=hvgs)
sce <- runUMAP(sce, dimred="PCA")

mat <- reducedDim(sce, "PCA")
dim(mat)

### library(bluster)
hclust.out <- clusterRows(mat, HclustParam())
plotUMAP(sce, colour_by=I(hclust.out))

hp2 <- HclustParam(method="ward.D2", cut.dynamic=TRUE)
hclust.out <- clusterRows(mat, hp2)
plotUMAP(sce, colour_by=I(hclust.out))

### Affinity propagation    library("apcluster")
set.seed(1000)
sub <- sce[,sample(ncol(sce), 200)]
ap.out <- clusterRows(reducedDim(sub), AffinityParam())
plotUMAP(sub, colour_by=I(ap.out))

set.seed(1000)
ap.out <- clusterRows(reducedDim(sub), AffinityParam(q=-2))
plotUMAP(sub, colour_by=I(ap.out))

### Centroid-based K-Means
set.seed(100)
kmeans.out <- clusterRows(mat, KmeansParam(10))
plotUMAP(sce, colour_by=I(kmeans.out))

kp <- KmeansParam(sqrt)
kp
set.seed(100)
kmeans.out <- clusterRows(mat, kp)
plotUMAP(sce, colour_by=I(kmeans.out))

set.seed(100)
mbkmeans.out <- clusterRows(mat, bluster::MbkmeansParam(20))
plotUMAP(sce, colour_by=I(mbkmeans.out))

### SOM kohonen
set.seed(1000)
som.out <- clusterRows(mat, SomParam(20))
plotUMAP(sce, colour_by=I(som.out))

set.seed(1000)
som.out <- clusterRows(mat, SomParam(100), full=TRUE)

par(mfrow=c(1,2))
plot(som.out$objects$som, "counts")
grid <- som.out$objects$som$grid$pts
text(grid[,1], grid[,2], seq_len(nrow(grid)))

### graph-based igraph
set.seed(101) # just in case there are ties.
graph.out <- clusterRows(mat, NNGraphParam(k=10))
plotUMAP(sce, colour_by=I(graph.out))

set.seed(101) # just in case there are ties.
np <- NNGraphParam(k=20, cluster.fun="louvain")
graph.out <- clusterRows(mat, np)
plotUMAP(sce, colour_by=I(graph.out))

### Density-based DBScan
dbscan.out <- clusterRows(mat, DbscanParam())
plotUMAP(sce, colour_by=I(dbscan.out))

summary(is.na(dbscan.out))

dbscan.out <- clusterRows(mat, DbscanParam(core.prop=0.1))
summary(is.na(dbscan.out))

### TWo-phase (multipleclustering methods overlapped for large datasets)
set.seed(100) # for the k-means
two.out <- clusterRows(mat, TwoStepParam())
plotUMAP(sce, colour_by=I(two.out))
 
twop <- TwoStepParam(second=NNGraphParam(k=5))
set.seed(100) # for the k-means
two.out <- clusterRows(mat, TwoStepParam())
plotUMAP(sce, colour_by=I(two.out))

### Statistics
nn.out <- clusterRows(mat, NNGraphParam(), full=TRUE)
nn.out$objects$graph

table(nn.out$clusters)


### Bluster package use is to easily control different algorithms. 
# The r Biocpkg("bluster") package provides a flexible and 
# extensible framework for clustering in Bioconductor packages/workflows. 
# At its core is the clusterRows() generic that controls dispatch to different 
# clustering algorithms. We will demonstrate on some single-cell RNA 
# sequencing data from the r Biocpkg("scRNAseq") package; 
# our aim is to cluster cells into cell populations based on their PC coordinates.

# clusterRows() enables users or developers to easily switch between clustering 
# algorithms by changing a single argument. Indeed, by passing the BlusterParam 
# object across functions, we can ensure that the same algorithm is used through 
# a workflow. It is also helpful for package functions as it provides diverse 
# functionality without compromising a clean function signature. However, the 
# true power of this framework lies in its extensibility. Anyone can write a 
# clusterRows() method for a new clustering algorithm with an associated 
# BlusterParam subclass, and that procedure is immediately compatible with any 
# workflow or function that was already using clusterRows().

### Bibliography
# https://rdrr.io/github/LTLA/bluster/f/vignettes/clusterRows.Rmd
# https://github.com/LTLA/bluster
# https://www.bioconductor.org/packages/release/bioc/manuals/bluster/man/bluster.pdf
