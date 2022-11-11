# title: "TSNE"
# author: "Casper"
# date: "19-10-2022"

### Load required packages
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(Rtsne)
  library(dgeAnalysis)
})

### TSNE
plot_data <- tsne_data(normDge)

### Scatterplot TSNE Celltype
scatter_plot(
  df = plot_data,
  x = "V1",
  y = "V2",
  group = "Celltype", 
  size = 5,
  title = "tSNE",
  xlab = "tSNE 1",
  ylab = "tSNE 2",
)

### Scatterplot TSNE Week
scatter_plot(
  df = plot_data,
  x = "V1",
  y = "V2",
  group = "Week", 
  size = 5,
  title = "tSNE",
  xlab = "tSNE 1",
  ylab = "tSNE 2",
)

### SessionInfo
sessionInfo()