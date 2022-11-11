library(tidyverse)
library(cluster)
library(factoextra)

df <- USArrests
df <- na.omit(df)
df <- scale(df)
head(df)

set.seed(123)

wss <- function(k) {
  kmeans(df, k, nstart=10)$tot.withinss
}

k.values <- 1:15

wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
                   type="b", pch=19, frame=FALSE,
                   xlab="Number of clusters K",
                   ylab="Total within-clusters sum of squares")

fviz_nbclust(df, kmeans, method="wss")

k4 <- kmeans(df, centers=4, nstart=25)
fviz_cluster(k4, df)
