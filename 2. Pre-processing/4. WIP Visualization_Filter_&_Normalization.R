# title: "Visualization Filter & Normalization"
# author: "Casper Kant"
# date: "19-10-2022"

### Load required packages


dge <- dge[rowSums(abs(dge$counts)) > 1, ]

### Gene count distribution plot before filtering
# tempDge <- dge
# tempDge$counts <- cpm(dge, log = TRUE, prior.count = 1)
# countDistributionLinePlot(tempDge)

### Filter genes
edgeR <- calcNormFactors(dge, method = "TMM")
cpm_perc = 25
cpm_value = 1
counts <- cpm(edgeR, log = TRUE, prior.count = 1)
selectedFeatures <- rownames(edgeR)[apply(counts, 1, function(v)
  sum(v >= cpm_value)) >= (cpm_perc / 100) * ncol(counts)]
highExprDge <- dge[selectedFeatures, , keep.lib.sizes = FALSE]

# Implement way to filter out household genes that have a constant expression (within SD) across all samples.

### Normalize
normDge <- calcNormFactors(highExprDge, method = "TMM")
normDge$counts <- cpm(normDge, log = TRUE, prior.count = 1)

### Gene count distribution plot after filtering
# tempDge <- normDge
# tempDge$counts <- cpm(normDge, log = TRUE, prior.count = 1)
# countDistributionLinePlot(tempDge)