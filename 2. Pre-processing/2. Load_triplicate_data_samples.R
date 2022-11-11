# title: "Load triplicate data samples"
# author: "Casper Kant"
# date: "19-10-2022"

### Load required packages
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
  library(dgeAnalysis)
})

### Read data
data_samples <- read.table(
  file = "C:/LUMC_Internship/Data/trisamples_KS2.csv",
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

### Take average of triplicate data
subset_data <- seq(from=1, to=60, by=3)
data_tricounts <- matrix(nrow= nrow(data_counts), ncol = ncol(data_counts)/3)
k = 1
for (i in subset_data) {
  data_tricounts[,k]  <- rowMeans(data_counts[c(i:(i+2))], na.rm = FALSE, dims = 1)
  k <- k+1
}

### Rename rows & columns
row.names(data_tricounts) <- as.character(rownames(data_counts))
colnames(data_tricounts) <- c("MSC0", "pLV0", "fl0", "tr0", "mut0",
                              "MSC1", "pLV1", "fl1", "tr1", "mut1",
                              "MSC2", "pLV2", "fl2", "tr2", "mut2",
                              "MSC3", "pLV3", "fl3", "tr3", "mut3")

### Read files as SE
se <- readCountsFromTable(data_tricounts[!grepl('^__', rownames(data_tricounts)), ], data_samples)
se <- addSamplesFromTableToSE(se, data_samples)

### Create DGE object
dge <- DGEList(counts = assay(se),
               samples = colData(se),
               genes = rowData(se))
row.names(dge$genes) <- row.names(dge$counts)

### Clean up Environment
rm(i, k, subset_data)

### Check loaded data
dim(dge$counts)  # 61541 20
dim(dge$samples) # 20 5
# View(dge$counts)
# View(dge$samples)

### SessionInfo
sessionInfo()