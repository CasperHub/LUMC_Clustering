### Import required packages  
suppressPackageStartupMessages({
  if (!require("kohonen")) install.packages("kohonen")
  devtools::install_github("LUMC/dgeAnalysis")
  library("dgeAnalysis")
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


### Plot PCA
plot_data <- pca_data(normDge) #plot with 9 PC's
pca_data <- plot_data[c(2:4, 14:16, 38:40)]
normDge$counts <- cpm(normDge, log = TRUE, prior.count = 1)
rm(dge, counts, cpm_perc, cpm_value, edgeR, highExprDge, selectedFeatures, se, plot_data)

### Kohonen stuff
som.data = som(scale(pca_data), grid = somgrid(xdim = 7, ydim = 7, "hexagonal"))
plot(som.data, type = "mapping")


### Kohnen plot help guide
data(wines)
set.seed(7)

kohmap <- xyf(scale(wines), vintages,
              grid = somgrid(5, 5, "hexagonal"), rlen=100)
plot(kohmap, type="changes")
counts <- plot(kohmap, type="counts", shape = "straight")

## show both sets of codebook vectors in the map
par(mfrow = c(1,2))
plot(kohmap, type="codes", main = c("Codes X", "Codes Y"))

par(mfrow = c(1,1))
similarities <- plot(kohmap, type="quality", palette.name = terrain.colors)
plot(kohmap, type="mapping",
     labels = as.integer(vintages), col = as.integer(vintages),
     main = "mapping plot")

## add background colors to units according to their predicted class labels
xyfpredictions <- classmat2classvec(getCodes(kohmap, 2))
bgcols <- c("gray", "pink", "lightgreen")
plot(kohmap, type="mapping", col = as.integer(vintages),
     pchs = as.integer(vintages), bgcol = bgcols[as.integer(xyfpredictions)],
     main = "another mapping plot", shape = "straight", border = NA)

## Show 'component planes'
set.seed(7)
sommap <- som(scale(wines), grid = somgrid(6, 4, "hexagonal"))
plot(sommap, type = "property", property = getCodes(sommap, 1)[,1],
     main = colnames(getCodes(sommap, 1))[1])

## Show the U matrix
Umat <- plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")

## use hierarchical clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(sommap, "codes")), 5)
add.cluster.boundaries(sommap, som.hc)

## and the same for rectangular maps
set.seed(7)
sommap <- som(scale(wines),grid = somgrid(6, 4, "rectangular"))
plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")

## use hierarchical clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(sommap, "codes")), 5)
add.cluster.boundaries(sommap, som.hc)

########################

data_train <- data[,c(2,4,5,8)]
 
data_train_matrix <- as.matrix(scale(data_train))

som_grid <- somgrid(xdim = 20, ydim=20, topo="hexagonal")

som_model <- som(data_train_matrix, 
                 grid=som_grid, 
                 rlen=100, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE,
                 n.hood="circular" )

### Practise code
library(scRNAseq)
sce <- ZeiselBrainData()

# Trusting the authors' quality control, and going straight to normalization.
library(scuttle)
sce <- logNormCounts(sce)

# Feature selection based on highly variable genes.
library(scran)
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n=1000)

# Dimensionality reduction for work (PCA) and pleasure (t-SNE).
set.seed(1000)
library(scater)
sce <- runPCA(sce, ncomponents=20, subset_row=hvgs)
sce <- runUMAP(sce, dimred="PCA")

mat <- reducedDim(sce, "PCA")
dim(mat)

set.seed(1000)
som.out <- clusterRows(mat, SomParam(20))
plotUMAP(sce, colour_by=I(som.out))

set.seed(1000)
som.out <- clusterRows(mat, SomParam(100), full=TRUE)

par(mfrow=c(1,2))
plot(som.out$objects$som, "counts")
grid <- som.out$objects$som$grid$pts
text(grid[,1], grid[,2], seq_len(nrow(grid)))

### Bibliography
# https://cran.r-project.org/web/packages/kohonen/kohonen.pdf

data(wines)

## som
som.wines <- som(scale(wines), grid = somgrid(5, 5, "hexagonal"))

## xyf
xyf.wines <- xyf(scale(wines), vintages, grid = somgrid(5, 5, "hexagonal"))

## supersom example
data(yeast)
yeast.supersom <- supersom(yeast, somgrid(6, 6, "hexagonal"),
                           whatmap = c("alpha", "cdc15", "cdc28", "elu"),
                           maxNA.fraction = .5)

plot(yeast.supersom, "changes")

obj.classes <- as.integer(yeast$class)
colors <- c("yellow", "green", "blue", "red", "orange")
plot(yeast.supersom, type = "mapping", col = colors[obj.classes],
     pch = obj.classes, main = "yeast data")



iris.sc = scale(iris[, 1:4])
iris.grid = somgrid(xdim = 5, ydim=5, topo="hexagonal")

# build model
iris.som = som(iris.sc, grid=iris.grid, rlen=100, alpha=c(0.05,0.01))

