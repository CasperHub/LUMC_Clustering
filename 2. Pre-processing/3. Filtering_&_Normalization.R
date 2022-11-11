# title: "Filtering & Normalization"
# author: "Casper Kant"
# date: "19-10-2022"

### Load required packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
  library(dgeAnalysis)
})

### Filter genes
dge <- dge[rowSums(abs(dge$counts)) > 1, ]
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

### Check loaded data
dim(normDge$counts) #12760 60
summary(normDge$counts)
# View(normDge$counts)

### SessionInfo
sessionInfo()