######################## IDEAS DISCUSSED ########################
#' - prereq for auc_exp4.2.R
#' - shows the validity of a hierarchical clustering-based classifier
#' - high accuracy on our limited training set, with accuracy improving as
#'   we cluster on the PC's and remove sparse drug entries 
#################################################################

source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")

#let's do clustering on two cancers w two groups to see something
#eventually scale this process to a plot with multiple cancer pairings 

sort(table(cancer_freqs), decreasing = TRUE)
cancers = c("BREAST", "CENTRAL_NERVOUS_SYSTEM")

cancer.data = block_dat(cancers, auc)

library(zoo)
cancer.data = na.aggregate(cancer.data)

#get samples as auc scores 
cancer.data = t(cancer.data)

### PCA
cancer.cov = cov(cancer.data)
cancer.eigen = eigen(cancer.cov)
scores = cancer.data%*%cancer.eigen$vectors

#clustering on un-projected
cancer.dist = dist(cancer.data, method="euclidean")
cancer.hclust = hclust(cancer.dist, method = "complete")

x11()
plot(cancer.hclust)
dev.off()

cluster.ec <- cutree(cancer.hclust, k=2) 

#now get `real` labels and colour maps
labels.real = str_util(rownames(cancer.data))
reference_map = ifelse(labels.real == "BREAST", 'red', 'blue')
cluster_map = ifelse(cluster.ec==1,'red','blue')

x11()
par(mfrow=c(1,2))
plot(scores[,1], scores[,2], col=cluster_map, pch=16, asp=1, main="hierarchical clustering", lwd=2)
plot(scores[,1], scores[,2], col=reference_map, pch=16, asp=1, main="real labels", lwd=2)
dev.off()

#quantify the results
metric = data.frame(cancer = labels.real, label = cluster.ec)
(table(metric)[1,1]+table(metric)[2,2])/sum(table(metric))


#clustering on pc's

#let's investigate the principle components some more 
x11()
par(mfrow=c(1,2))
plot(cancer.eigen$values, xlab="pc", ylab="variance")
plot(cumsum(cancer.eigen$values)/sum(cancer.eigen$values), xlab="pc", ylab="explained variance")
abline(h=0.85, col='red', lty=2)
dev.off()

x11()
par(mfrow=c(3,1))
barplot(cancer.eigen$vectors[,1])
barplot(cancer.eigen$vectors[,2])
barplot(cancer.eigen$vectors[,3])
dev.off()

#select number of principle components corresponding to 90% of explained variance 
n.comp = min(which(cumsum(cancer.eigen$values)/sum(cancer.eigen$values)>=0.9))
n.comp

#now cluster again 
cancer.reduced = scores[,1:n.comp]
reduced.dist = dist(cancer.reduced, method="euclidean")
reduced.hclust = hclust(reduced.dist, method = "complete")

x11()
plot(reduced.hclust)
dev.off()

reduced.ec <- cutree(reduced.hclust, k=2) 

#plot
cluster_map.reduced = ifelse(reduced.ec==1,'red','blue')

x11()
par(mfrow=c(1,2))
plot(scores[,1], scores[,2], col=cluster_map.reduced, pch=16, asp=1, main="hierarchical clustering", lwd=2)
plot(scores[,1], scores[,2], col=reference_map, pch=16, asp=1, main="real labels", lwd=2)
dev.off()

metric.reduced = data.frame(cancer = labels.real, label = reduced.ec)
table(metric.reduced)
(table(metric.reduced)[1,1]+table(metric.reduced)[2,2])/sum(table(metric.reduced))


#WE NOW REDUCE THE TOTAL NA COUNT ON A DRUGWISE BASIS 

#Now let's do the same type of thing with fewer treatments 
cancer2.data = block_dat(cancers, auc)
na_rows = sort(rowSums(is.na(cancer2.data)), decreasing=TRUE)

#we can choose to reduce num rows such that ratio of nas to sample size is larger than 50%
num_samples = dim(cancer2.data)[2]
drugs_kept = names(which(na_rows/num_samples <0.5))

#reduce 
cancer2.data = cancer2.data[drugs_kept, ]

#go through the clustering once again 
cancer2.data = na.aggregate(cancer2.data)
cancer2.data = t(cancer2.data)

### PCA
cancer2.cov = cov(cancer2.data)
cancer2.eigen = eigen(cancer2.cov)

#visualization
x11()
par(mfrow=c(1,2))
plot(cancer2.eigen$values, xlab="pc", ylab="variance")
plot(cumsum(cancer2.eigen$values)/sum(cancer2.eigen$values), xlab="pc", ylab="explained variance")
abline(h=0.85, col='red', lty=2)
dev.off()

#select some principal components 
n2.comp = min(which(cumsum(cancer2.eigen$values)/sum(cancer2.eigen$values)>=0.9))
n2.comp

scores2 = cancer2.data%*%cancer2.eigen$vectors
cancer2.reduced = scores2[,1:n2.comp]

#clustering on un-projected
cancer2.dist = dist(cancer2.data, method="euclidean")
cancer2.hclust = hclust(cancer2.dist, method = "complete")

x11()
plot(cancer2.hclust)
dev.off()

cluster2.ec <- cutree(cancer2.hclust, k=2) 

#now get `real` labels and colour maps
labels2.real = str_util(rownames(cancer2.data))
reference_map2 = ifelse(labels2.real == "BREAST", 'red', 'blue')
cluster_map2 = ifelse(cluster2.ec==1,'red','blue')

x11()
par(mfrow=c(1,2))
plot(scores2[,1], scores2[,2], col=cluster_map2, pch=16, asp=1, main="hierarchical clustering", lwd=2)
plot(scores2[,1], scores2[,2], col=reference_map2, pch=16, asp=1, main="real labels", lwd=2)
dev.off()

metric2 = data.frame(cancer = labels2.real, label = cluster2.ec)
table(metric2)
(table(metric2)[1,1]+table(metric2)[2,2])/sum(table(metric2))

#now let's try with the first few pc's 
reduced2.dist = dist(cancer2.reduced, method="euclidean")
reduced2.hclust = hclust(reduced2.dist, method = "complete")

x11()
plot(reduced2.hclust)
dev.off()

reduced2.ec <- cutree(reduced2.hclust, k=2) 

#plot
cluster_map.reduced2 = ifelse(reduced2.ec==1,'red','blue')

x11()
par(mfrow=c(1,2))
plot(scores2[,1], scores2[,2], col=cluster_map.reduced2, pch=16, asp=1, main="hierarchical clustering", lwd=2)
plot(scores2[,1], scores2[,2], col=reference_map2, pch=16, asp=1, main="real labels", lwd=2)
dev.off()

metric.reduced2 = data.frame(cancer = labels2.real, label = reduced2.ec)
table(metric.reduced2)
(table(metric.reduced2)[1,1]+table(metric.reduced2)[2,2])/sum(table(metric.reduced2))

#okay this is a bit more accurate 

#color map for the misclassified only?
difference_map = rep('grey', dim(scores2)[1])
difference_map[which(cluster_map.reduced2!=reference_map2)] = 'black'

#can check accuracy this way too 
p = length(which(cluster_map.reduced2!=reference_map2))/length(reference_map2)
max(p,1-p)

x11()
plot(scores2[,1], scores2[,2], col=difference_map, pch=16, asp=1, main="hierarchical clustering", lwd=2)
dev.off()

#i've noticed increasing accuracy from all data no pc -> all data pc -> reduced data no pc -> ... 
