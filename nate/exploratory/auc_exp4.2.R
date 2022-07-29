source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")
library(zoo)

#let's do clustering on two cancers w two groups to see something
#eventually scale this process to a plot with multiple cancer pairings 

cancers = sort(table(cancer_freqs), decreasing = TRUE)[1:10]
cancers = names(cancers)
cancer_pairs = combn(cancers, 2)

accuracies = c(NULL)
pairs = c(NULL)

for(i in 1:dim(cancer_pairs)[2]){
  pair = cancer_pairs[,i]
  
  cancer.data = block_dat(pair, auc)
  na_rows = sort(rowSums(is.na(cancer.data)), decreasing=TRUE)
  
  #we can choose to reduce num rows such that ratio of nas to sample size is larger than 50%
  num_samples = dim(cancer.data)[2]
  drugs_kept = names(which(na_rows/num_samples <0.5))
  
  #reduce 
  cancer.data = cancer.data[drugs_kept, ]
  
  #go through the clustering once again 
  cancer.data = na.aggregate(cancer.data)
  cancer.data = t(cancer.data)
  
  ### PCA
  cancer.cov = cov(cancer.data)
  cancer.eigen = eigen(cancer.cov)
  
  #visualization
  #' x11()
  #' par(mfrow=c(1,2))
  #' plot(cancer.eigen$values, xlab="pc", ylab="variance")
  #' plot(cumsum(cancer.eigen$values)/sum(cancer.eigen$values), xlab="pc", ylab="explained variance")
  #' abline(h=0.85, col='red', lty=2)
  #' dev.off()
  
  #select some principal components 
  n.comp = min(which(cumsum(cancer.eigen$values)/sum(cancer.eigen$values)>=0.9))
  n.comp
  
  scores = cancer.data%*%cancer.eigen$vectors
  cancer.reduced = scores[,1:n.comp]
  
  #now get `real` labels and colour maps
  labels.real = str_util(rownames(cancer.data))
  reference_map = ifelse(labels.real == "BREAST", 'red', 'blue')
  
  #now let's try with the first few pc's 
  reduced.dist = dist(cancer.reduced, method="euclidean")
  reduced.hclust = hclust(reduced.dist, method = "complete")
  
  reduced.ec <- cutree(reduced.hclust, k=2) 
  #plot
  cluster_map.reduced = ifelse(reduced.ec==1,'red','blue')
  
  #' x11()
  #' par(mfrow=c(1,2))
  #' plot(scores[,1], scores[,2], col=cluster_map.reduced, pch=16, asp=1, main="hierarchical clustering", lwd=2)
  #' plot(scores[,1], scores[,2], col=reference_map, pch=16, asp=1, main="real labels", lwd=2)
  #' dev.off()
  
  metric.reduced = data.frame(cancer = labels.real, label = reduced.ec)
  table(metric.reduced)
  accuracy = (table(metric.reduced)[1,1]+table(metric.reduced)[2,2])/sum(table(metric.reduced))
  
  accuracies = c(accuracies, accuracy)
  pairs = c(pairs, paste(pair[1], " ", pair[2]))
}

names(accuracies)<-pairs





