# title: "DBSCAN"
# author: "Casper"
# date: "24-10-2022"

### Packages
library(fpc)
library(dbscan)
library(factoextra)

### Data import
data("iris")
new <- iris[, -5]

### Obtaining eps value
kNNdistplot(new, k=3)
abline(h=.45, lty=2)

### Clustering
set.seed(123)
f <- fpc::dbscan(new, eps=0.45, MinPts = 4)
d <- dbscan::dbscan(new, eps = 0.45, minPts = 4)

### Visualization
fviz_cluster(f, new, geom = "points")  
fviz_cluster(d, new, geom = "points")
