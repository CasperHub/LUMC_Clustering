# title: "Connectivity-based Hierarchical"
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

# Fitting Hierarchical clustering Model to training dataset
set.seed(123)
Hierar_cl <- hclust(D, method = "average")

# Plotting dendrogram
Hierar_cl$labels <- colnames(normDge$counts)
plot(Hierar_cl)

# Choosing no. of clusters cutting tree by height
# abline(h = 40, col = "green")

# Cutting tree by no. of clusters
cluster <- cutree(Hierar_cl, k = 10 )
cluster

table(cluster)
rect.hclust(Hierar_cl, k = 10, border = "green")

### Obtain quantitative score silhouette
results <- matrix(nrow= 9, ncol = 2)
for (i in 2:10) {
  set.seed(123)
  cluster <- cutree(Hierar_cl, k = i )
  results[i-1] <- mean(data.matrix(silhouette(cluster, D)[c(121:180)]))
  plot(silhouette(cluster, D), col=1:i) 
}

colnames(results) <- c("Scores", "Clusters")
results[10:18] <- c(2:10)
results <- data.frame(results)
results[rev(order(results$Scores)),] 
results[which.max(results$Scores),]

### Scoring on single cluster result
# cluster <- data.matrix(cluster)
# plot(silhouette(cluster, D), col=1:6)
# score <- silhouette(cluster, D)
# score <- data.frame(score)
# sum(score$sil_width) / length(score$sil_width)

### Plotting formed clusters on PCA with interactivity
plot_data <- pca_data(normDge)
set.seed(123)
cluster <- cutree(Hierar_cl, k = 10 )
clusters <- data.frame(clusters = as.character(cluster))
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
  layout(title="Hierarchical PCA clustering")
p1

### Save to HTML plot
saveWidget(p1, "Hierarchical_PCA_interactive.html", selfcontained=FALSE)

### Plotting formed clusters on TSNE with interactivity
plot_data <- tsne_data(normDge)
set.seed(123)
cluster <- cutree(Hierar_cl, k = 10 )
clusters <- data.frame(clusters = as.character(cluster))
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
  theme_bw() + #remove + to add legend
  theme(legend.position='none')

p3 <- ggplotly(p2, tooltip = "text") %>% 
  layout(title="Hierarchical tSNE clustering")
p3

### Save to HTML plot
saveWidget(p3, "Hierarchical_tSNE_interactive.html", selfcontained=FALSE)