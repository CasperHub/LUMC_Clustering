# title: "Full_PCA_Centroid_Kmeans_Silhouette"
# author: "Casper"
# date: "19-10-2022"

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

### K-means + decide K
set.seed(123)
k = 4
km = kmeans(pca_data, centers = k, nstart = 50)

### Plotting cluster centers
km$centers
km$centers[, c("PC1", "PC2")]

### Model Evaluation and visualization
plot(pca_data[c("PC1", "PC2")],
     col = km$cluster,
     main = "K-means with 5 clusters")

### cex is font size, pch is symbol
points(km$centers[, c("PC1", "PC2")],
       col = 1:5, pch = 8, cex = 3)

### Visualizing clusters
clusplot(pca_data[, c("PC1", "PC2")],
         km$cluster,
         lines = 0,
         shade = TRUE,
         color = TRUE,
         plotchar = TRUE,
         span = TRUE,
         main = paste("K-means clustering PCA"),
         xlab = 'PC1',
         ylab = 'PC2'
)

### Obtain quantitative score (sillhouette)
D <- daisy(pca_data)
results <- matrix(nrow= 9, ncol = 2)

for (i in 2:10) {
  set.seed(123)
  km = kmeans(pca_data, centers = i, nstart = 50)
  results[i-1] <- mean(data.matrix(silhouette(km$cluster, D)[c(121:180)]))
  plot(silhouette(km$cluster, D), col=1:i)
}

colnames(results) <- c("Scores", "Clusters")
results[10:18] <- c(2:10)
results <- data.frame(results)
results[rev(order(results$Scores)),] 
results[which.max(results$Scores),]

### Interactivity plot PCA cluster
plot_data <- pca_data(normDge)
set.seed(123)
km = kmeans(pca_data, centers = 5, nstart = 50)
clusters <- data.frame(clusters = as.character(km$cluster))
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
  layout(title="K-means PCA clustering")
p1

### Save to HTML plot
saveWidget(p1, "K-means_PCA_interactive.html", selfcontained=FALSE)

### Interactivity plot tSNE cluster
plot_data <- tsne_data(normDge)
set.seed(123)
km = kmeans(pca_data, centers = 5, nstart = 50)
clusters <- data.frame(clusters = as.character(km$cluster))
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
  layout(title="K-means tSNE clustering")
p3

### Save to HTML plot
saveWidget(p3, "K-means_tSNE_interactive.html", selfcontained=FALSE)