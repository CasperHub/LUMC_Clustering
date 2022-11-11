# norm.counts_KS2 <- read.delim("C:/LUMC_Internship/Data/norm.counts_KS2.txt")
# week1 <- read.delim("C:/LUMC_Internship/Data/week1.txt")
# week2 <- read.delim("C:/LUMC_Internship/Data/week2.txt")
# week3 <- read.delim("C:/LUMC_Internship/Data/week3.txt")

library(tidyverse)
library(cluster)
library(factoextra)

week0 <- read.delim("C:/LUMC_Internship/Data/week0.txt", row.names=1)
week0 <- na.omit(week0)
week0 <- scale(week0)
head(week0)

set.seed(123)

# wss <- function(k) {
#   kmeans(week0, k, nstart=10)$tot.withinss
# }
# 
# k.values <- 1:15
# 
# wss_values <- map_dbl(k.values, wss)
# 
# plot(k.values, wss_values,
#      type="b", pch=19, frame=FALSE,
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")

fviz_nbclust(week0, kmeans, method="wss")

# avg_sil <- function(k) {
#   km.res <- kmeans(week0, centers=k,nstart=25)
#   ss <- silhouette(km.res$cluster,dist(week0))
#   mean(ss[, 3])
# }
# 
# k.values <-2:15
# 
# avg_sil_values <- map_dbl(k.values, avg_sil)
# 
# plot(k.values, avg_sil_values, type="b", pch=19, frame=FALSE, xlab="Number of clusters K", ylab="Average Silhouettes")

fviz_nbclust(week0, kmean, method="silhouette")
# fviz_gap_stat(gap_stat)