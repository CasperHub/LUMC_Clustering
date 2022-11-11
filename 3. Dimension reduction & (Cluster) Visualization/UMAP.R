# title: "TSNE"
# author: "Casper"
# date: "19-10-2022"

### Load required packages
suppressPackageStartupMessages({
  library(umap)
  library(dplyr)
})

### UMAP
umap_results <- umap::umap(normDge$counts, intgroup=c("SUBTYPE"), returnData=TRUE)
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("PATIENT_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(clinical, by = "PATIENT_ID")
  ggplot(umap_plot_df,aes(x = X1,y = X2,
                        color = SUBTYPE,
                        shape = DSS_STATUS)) + geom_point()


# ### SessionInfo
# sessionInfo()