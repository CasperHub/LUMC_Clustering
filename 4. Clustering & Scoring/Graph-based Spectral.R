# title: "Graph-based Spectral"
# author: "Casper"
# date: "10-11-2022"

### Import required packages
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  library(rlang)
  library(limma)
  library(edgeR)
  library(cluster)
  library(ggplot2)
  library(factoextra)
  library(SummarizedExperiment)
  library(dgeAnalysis)
  library(Rtsne)
  library(plotly)
  library(htmlwidgets)
})

### Clustering scoring metrics pca
fviz_nbclust(pca_data, kmeans, method = 'wss')
fviz_nbclust(pca_data, kmeans, method = 'silhouette')

# Finding distance matrix
D <- daisy(pca_data) 
#distance_mat <- dist(pca_data, method = 'euclidean') //Alternative distance matrix option

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

# ### Visualizing Spectral clustering
# scatter_plot(
#   df = X_sc,
#   x = "X1",
#   y = "X2",
#   group = "cluster",
#   size = 4,
#   title = "Spectral",
# )

### Silhouette Scoring
# D <- daisy(pca_data)
# plot(silhouette(X_sc_kmeans$cluster, D), col=1:4)
# score <- silhouette(X_sc_kmeans$cluster, D)
# score <- data.frame(score)
# sum(score$sil_width) / length(score$sil_width)

### Obtain quantitative score silhouette
results <- matrix(nrow= 9, ncol = 2)
for (i in 2:10) {
  set.seed(123)
  X_sc_kmeans <- kmeans(X_sc, i)
  results[i-1] <- mean(data.matrix(silhouette(X_sc_kmeans$cluster, D)[c(121:180)]))
  plot(silhouette(X_sc_kmeans$cluster, D), col=1:i)
}

colnames(results) <- c("Scores", "Clusters")
results[10:18] <- c(2:10)
results <- data.frame(results)
results[rev(order(results$Scores)),] 
results[which.max(results$Scores),]


### Plotting formed clusters on PCA with interactivity
plot_data <- pca_data(normDge)
set.seed(123)
X_sc_kmeans <- kmeans(X_sc, 3)
clusters <- data.frame(clusters = as.character(X_sc_kmeans$cluster))
plot_data <- cbind(plot_data, clusters)

p <- ggplot(plot_data,
            aes(x=PC1,
                y=PC2,
                #alpha=Week,
                color=clusters,
                shape=Celltype,
                text=paste0("Sample:", sample,
                            "\nCelltype: ", Celltype,
                            "\nCluster: ", clusters,
                            "\nWeek: ", Week))) +
  geom_point(size=4) +
  theme_bw() + #remove + to add legend
  theme(legend.position='none')

p1 <- ggplotly(p, tooltip = "text") %>%
  layout(title="Spectral PCA clustering")
p1

### Save to HTML plot
saveWidget(p1, "Spectral_PCA_interactive.html", selfcontained=FALSE)

### Plotting formed clusters on TSNE with interactivity
plot_data <- tsne_data(normDge)
set.seed(123)
X_sc_kmeans <- kmeans(X_sc, 3)
clusters <- data.frame(clusters = as.character(X_sc_kmeans$cluster))
plot_data <- cbind(plot_data, clusters)

p2 <- ggplot(plot_data,
             aes(x=V1,
                 y=V2,
                 color=clusters,
                 shape=Celltype,
                 text=paste0("Sample:", sample,
                             "\nCelltype: ", Celltype,
                             "\nCluster: ", clusters,
                             "\nWeek: ", Week))) +
  geom_point(size=4) +
  theme_bw() +
  theme(legend.position='none')

p3 <- ggplotly(p2, tooltip = "text") %>%
  layout(title="Spectral tSNE clustering")
p3

### Save to HTML plot
saveWidget(p3, "Spectral_tSNE_interactive.html", selfcontained=FALSE)