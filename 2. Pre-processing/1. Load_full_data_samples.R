# title: "Load full data samples"
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

### Check loaded data
dim(dge$counts)  # 61541 60
dim(dge$samples) # 60 6
# View(dge$counts)
# View(dge$samples)

### SessionInfo
sessionInfo()