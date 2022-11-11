# title: "KS2_Filter"
# author: "Casper"
# date: "30-09-2022"

## Load required libraries
if (!require("edgeR")) install.packages("edgeR")
if (!require("limma")) install.packages("limma")

### Read data
data_location = "C:/LUMC_Internship/Data/norm.counts_KS2.txt"
data <- read.delim(data_location, row.names=1)
data[61] <- NULL

### Create DGE object
dge <- DGEList(data)
dge <- dge[rowSums(abs(dge$counts)) > 1,]

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
normDge <- cpm(normDge, log = TRUE, prior.count = 1)

### PCA after filter
df = t(normDge)
pc <- prcomp(df, center=TRUE, scale=TRUE)
summary(pc)
