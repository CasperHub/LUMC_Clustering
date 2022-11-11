# title: "Full_PCA_Connectivity_Hierarchical_Silhouette"
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
  if (!require("ClusterR")) install.packages("ClusterR")
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("factoextra")) install.packages("factoextra")
  if (!require("SummarizedExperiment")) install.packages("SummarizedExperiment")
  library("dgeAnalysis")
  library(dplyr)
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
counts <- edgeR::cpm(edgeR, log = TRUE, prior.count = 1)
selectedFeatures <- rownames(edgeR)[apply(counts, 1, function(v)
  sum(v >= cpm_value)) >= (cpm_perc / 100) * ncol(counts)]
highExprDge <- dge[selectedFeatures, , keep.lib.sizes = FALSE]

### Normalize
normDge <- calcNormFactors(highExprDge, method = "TMM")
normDge$counts <- edgeR::cpm(normDge, log = TRUE, prior.count = 1)

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

### Filter PC's & rename columns
pca_data <- plot_data[c(2:4, 14:16, 38:40)]
rownames(pca_data) <- plot_data$sample

### Clustering scoring metrics pca
fviz_nbclust(pca_data, kmeans, method = 'wss')
fviz_nbclust(pca_data, kmeans, method = 'silhouette')

# Finding distance matrix
# D <- daisy(pca_data) //Alternative distance matrix option
distance_mat <- dist(pca_data, method = 'euclidean')

# Fitting Hierarchical clustering Model to training dataset
set.seed(123)
Hierar_cl <- hclust(distance_mat, method = "average")

# Plotting dendrogram
Hierar_cl$labels <- colnames(normDge$counts)
plot(Hierar_cl)

# Choosing no. of clusters cutting tree by height
abline(h = 40, col = "green")

# Cutting tree by no. of clusters
fit <- cutree(Hierar_cl, k = 6 )
fit

table(fit)
rect.hclust(Hierar_cl, k = 5, border = "green")

### Show plot for different K
# plot(Hierar_cl)
# rect.hclust(Hierar_cl, k = 4, border = "green")
# plot(Hierar_cl)
# rect.hclust(Hierar_cl k = 8, border = "green")

### Obtain quantitavtie score silhouette
# fit <- data.matrix(fit)
# plot(silhouette(fit, distance_mat), col=1:6)
# score <- silhouette(fit, distance_mat)
# score <- data.frame(score)
# sum(score$sil_width) / length(score$sil_width)

results <- matrix(nrow= 9, ncol = 2)
for (i in 2:10) {
  set.seed(123)
  fit <- cutree(Hierar_cl, k = i )
  results[i-1] <- mean(matrix(matrix(silhouette(fit, distance_mat))[c(121:180)]))
  plot(silhouette(fit, distance_mat), col=1:i)
}

colnames(results) <- c("Scores", "Clusters")
results[10:18] <- c(2:10)
results <- data.frame(results)
results[rev(order(results$Scores)),] 
results[which.max(results$Scores),]