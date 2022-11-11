# title: "PCA"
# author: "Casper"
# date: "19-10-2022"

### Load required packages
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(dgeAnalysis)
})

### PCA
plot_data <- pca_data(normDge)

### Scatterplot PCA Celltype
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

### Scatterplot PCA Week
scatter_plot(
  df = plot_data,
  x = "PC1",
  y = "PC2",
  group = "Week",
  size = 5,
  title = "PCA",
  xlab = paste0("PC1", " (", plot_data$percent[13], "%)"),
  ylab = paste0("PC2", " (", plot_data$percent[14], "%)"),
)

### Get highest PCs from PCA
# pc <- matrix(plot_data$percent)
# pca_data <- plot_data[c(2:4, 14:16, 38:40)]

### SessionInfo
sessionInfo()