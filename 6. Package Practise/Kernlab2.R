library(kernlab)
## Create toy data
x <- rbind(matrix(rnorm(10),,2),matrix(rnorm(10,mean=3),,2))
y <- matrix(c(rep(1,5),rep(-1,5)))
### Use as.kernelMatrix to label the cov. matrix as a kernel matrix
### which is eq. to using a linear kernel
K <- as.kernelMatrix(crossprod(t(x)))
K
svp2 <- ksvm(K, y, type="C-svc")
svp2 