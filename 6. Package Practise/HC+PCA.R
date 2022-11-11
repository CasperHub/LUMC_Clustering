library(tidyverse)
library(ggplot2)
library(dplyr)
library(cluster)
library(factoextra)

df <- read.delim("C:/LUMC_Internship/Data/norm.counts_KS2.txt", row.names=1)
df[61] <- NULL
row_sub = apply(df, 1, function(row) all(row !=0 ))
df = df[row_sub,]
df = t(df)

pc <- prcomp(df, center=TRUE, scale=TRUE)
summary(pc)
df_transform = as.data.frame(-pc$x[,1:2])

hc <- df_transform %>%   #Get Data
      dist %>%  #Distance/Dissimilarity Matrix Euc
      hclust    #Compute Hc's

plot(hc)
