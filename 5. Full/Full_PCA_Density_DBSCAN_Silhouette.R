# title: "Full_PCA_Density_DBSCAN_Silhouette"
# author: "Casper"
# date: "21-10-2022"

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
  library(fpc)
  library(dbscan)
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
pca_data <- plot_data[c(2:4, 14:16, 38:40)]

# ### Clustering scoring metrics pca
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

### Obtain quantitavtie score silhouette
D <- daisy(pca_data)
plot(silhouette(d$cluster, D), col=1:5)

score <- silhouette(d$cluster, D)
score <- data.frame(score)
sum(score$sil_width) / length(score$sil_width)
