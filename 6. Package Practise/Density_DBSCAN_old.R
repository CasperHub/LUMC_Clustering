install.packages("dbscan")
install.packages("fpc")

library(tidyverse) 
library(cluster)
library(fpc)
library(dbscan)
library(factoextra)
library(ggplot2)
library(dplyr)

# df <- read.delim("C:/LUMC_Internship/Data/ks2_filtered.txt", row.names=1, sep=',')

df <- read.delim("C:/LUMC_Internship/Data/norm.counts_KS2.txt", row.names=1)
df[61] <- NULL
df = df[c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58)]
row_sub = apply(df, 1, function(row) all(row !=0 ))
df = df[row_sub,]
df = t(df)

pc <- prcomp(df, center=TRUE, scale=TRUE)
summary(pc)
new_df = as.data.frame(-pc$x[,1:2])

# scaling the dataset
new_df = scale(new_df)
new_df %>% head()
# 
# # to plot the eps values
# eps_plot = kNNdistplot(new_df, k=3)
# 
# # to draw an optimum line
# eps_plot %>% abline(h = 0.45, lty = 2)

set.seed(50)

# creation of an object km which store the output of the function kmeans
d <- dbscan::dbscan(new_df, eps = 0.45, MinPts =  2)
d

# cluster visualisation
fviz_cluster(d, new_df, geom = "point")

