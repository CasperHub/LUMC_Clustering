# title: "Density-based DBSCAN"
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
  library(fpc)
  library(dbscan)
  library(Rtsne)
  library(plotly)
  library(htmlwidgets)
})

### Clustering scoring metrics pca
fviz_nbclust(pca_data, kmeans, method = 'wss')
fviz_nbclust(pca_data, kmeans, method = 'silhouette')

### Checking optimal eps value
kNNdistplot(pca_data, k=5)
abline(h=46, lty=2)

### Clustering
set.seed(123)
# f <- fpc::dbscan(pca_data, eps=46, MinPts = 4)
d <- dbscan::dbscan(pca_data, eps = 46, minPts = 2)
fviz_cluster(d, pca_data, geom = "points")

### Obtain quantitative score silhouette
D <- daisy(pca_data)
plot(silhouette(d$cluster, D), col=1:5)

score <- data.frame(silhouette(d$cluster, D))
sum(score$sil_width) / length(score$sil_width)


### Plotting formed clusters on PCA with interactivity
plot_data <- pca_data(normDge)
set.seed(123)
d <- dbscan::dbscan(pca_data, eps = 46, minPts = 2)
clusters <- data.frame(clusters = as.character(d$cluster))
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
  layout(title="DBSCAN PCA clustering")
p1

### Save to HTML plot
saveWidget(p1, "DBSCAN_PCA_interactive.html", selfcontained=FALSE)

### Plotting formed clusters on TSNE with interactivity
plot_data <- tsne_data(normDge)
set.seed(123)
d <- dbscan::dbscan(pca_data, eps = 46, minPts = 2)
clusters <- data.frame(clusters = as.character(d$cluster))
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
  layout(title="DBSCAN tSNE clustering")
p3

### Save to HTML plot
saveWidget(p3, "DBSCAN_tSNE_interactive.html", selfcontained=FALSE)