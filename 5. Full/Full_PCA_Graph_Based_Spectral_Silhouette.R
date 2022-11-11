# title: "Full_PCA_Graph_Based_Spectral_Silhouette"
# author: "Casper"
# date: "19-10-2022"

### Import required packages
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!require("rlang")) install.packages("rlang")
  if (!require("limma")) install.packages("limma")
  if (!require("edgeR")) install.packages("edgeR")
  if (!require("cluster")) install.packages("cluster")
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("factoextra")) install.packages("factoextra")
  if (!require("SummarizedExperiment")) install.packages("SummarizedExperiment")
  library("dgeAnalysis")
  library(kernlab)
})

### Read data
data_samples <- read.table(
  file = "C:/LUMC_Internship/Data/samples_KS2.csv",
  header = TRUE,
  sep = ";",
  row.names = 1,
  check.names = FALSE
)

data_counts <- read.table(
  file = "C:/LUMC_Internship/Data/norm.counts_KS2.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
data_counts[61] <- NULL

### Read files as SE
se <- readCountsFromTable(data_counts[!grepl('^__', rownames(data_counts)), ], data_samples)
se <- addSamplesFromTableToSE(se, data_samples)

### Create DGE object
dge <- DGEList(counts = assay(se),
               samples = colData(se),
               genes = rowData(se))
row.names(dge$genes) <- row.names(dge$counts)
dge <- dge[rowSums(abs(dge$counts)) > 1, ]

### Filter genes
edgeR <- calcNormFactors(dge, method = "TMM")
cpm_perc = 25
cpm_value = 1
counts <- cpm(edgeR, log = TRUE, prior.count = 1)
selectedFeatures <- rownames(edgeR)[apply(counts, 1, function(v)
  sum(v >= cpm_value)) >= (cpm_perc / 100) * ncol(counts)]
highExprDge <- dge[selectedFeatures, , keep.lib.sizes = FALSE]

### Normalize
normDge <- calcNormFactors(highExprDge, method = "TMM")
normDge$counts <- cpm(normDge, log = TRUE, prior.count = 1)

### Plot PCA
plot_data <- pca_data(normDge) #plot with 9 PC's
scatter_plot(
  df = plot_data,
  x = "PC1",
  y = "PC2",
  group = "Celltype",
  size = 5,
  title = "PCA",
  xlab = paste0("PC1", " (", plot_data$percent[13], "%)"),
  ylab = paste0("PC2", " (", plot_data$percent[14], "%)"),
)

### Get highest PCs from PCA
pca_data <- plot_data[c(2:4, 14:16, 38:40)]

### Clustering scoring metrics pca
fviz_nbclust(pca_data, kmeans, method = 'wss')
fviz_nbclust(pca_data, kmeans, method = 'silhouette')

### Spectral clustering
spectral_clustering <- function(pca_data, # matrix of data points
                                nn = 10, # the k nearest neighbors to consider
                                n_eig = 2) # m number of eignenvectors to keep
{
  mutual_knn_graph <- function(pca_data, nn = 10)
  {
    # matrix of euclidean distances between data points in X
    D <- as.matrix( dist(pca_data) ) 
    
    # intialize the knn matrix
    knn_mat <- matrix(0,
                      nrow = nrow(pca_data),
                      ncol = nrow(pca_data))
    
    # find the 10 nearest neighbors for each point
    for (i in 1: nrow(pca_data)) {
      neighbor_index <- order(D[i,])[2:(nn + 1)]
      knn_mat[i,][neighbor_index] <- 1 
    }
    
    # Now we note that i,j are neighbors if K[i,j] = 1 or K[j,i] = 1 
    knn_mat <- knn_mat + t(knn_mat) # find mutual knn
    
    knn_mat[ knn_mat == 2 ] = 1
    
    return(knn_mat)
  }
  
  graph_laplacian <- function(W, normalized = TRUE)
  {
    stopifnot(nrow(W) == ncol(W)) 
    
    g = colSums(W) # degrees of vertices
    n = nrow(W)
    
    if(normalized)
    {
      D_half = diag(1 / sqrt(g) )
      return( diag(n) - D_half %*% W %*% D_half )
    }
    else
    {
      return( diag(g) - W )
    }
  }
  
  W = mutual_knn_graph(pca_data) # 1. matrix of similarities
  L = graph_laplacian(W) # 2. compute graph laplacian
  ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  n = nrow(L)
  return(ei$vectors[,(n - n_eig):(n - 1)]) # return the eigenvectors of the n_eig smallest eigenvalues
  
}

### do spectral clustering procedure
set.seed(123)
X_sc <- spectral_clustering(pca_data)

### run kmeans on the 4 eigenvectors
X_sc_kmeans <- kmeans(X_sc, 4)
X_sc <- data.frame(X_sc)

clusplot(X_sc[, c("X1", "X2")],
         X_sc_kmeans$cluster,
         lines = 0,
         shade = FALSE,
         color = FALSE,
         plotchar = TRUE,
         span = TRUE,
         main = paste("Spectral clustering PCA"),
         xlab = 'PC1',
         ylab = 'PC2'
)

### Alternative Spectral clustering kernlab
# sc <- specc(X, centers=5)
# plot(X, col=sc, pch=4)            # estimated classes (x)
# points(X, col=obj$classes, pch=5) # true classes 

X_sc <- cbind(cluster = matrix(X_sc_kmeans$cluster), X_sc)

# ## Visualizing Spectral clustering
# scatter_plot(
#   df = X_sc,
#   x = "X1",
#   y = "X2",
#   group = "cluster",
#   size = 4,
#   title = "Spectral",
# )

### Silhouette Scoring
# D <- daisy(pca_data) //Calculates own Dissimilarity matrix in spectral_clustering function. 
plot(silhouette(X_sc_kmeans$cluster, D), col=1:4)
score <- silhouette(X_sc_kmeans$cluster, D)
score <- data.frame(score)
sum(score$sil_width) / length(score$sil_width)
