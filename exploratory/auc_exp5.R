######################## IDEAS DISCUSSED ########################
#' - proper attempt at a pca on the drugs + clustering
#' - generate B/(B+W) curve to see optimal number
#' - plot our labelled drugs in the cluters on projected space 
#' - consider effect of pca 
#################################################################

#' NOTE: if we consider clustering on first n pc's we effectively say we only need to 
#' test drugs on specific loadings of cells!
#' 

source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")

#borrowing syntax from `auc_exp4.R`
cancers = cancer_types
cancer.data = block_dat(cancers, auc)
dim(cancer.data)

na_rows = sort(rowSums(is.na(cancer.data)), decreasing=TRUE)
num_cells = dim(cancer.data)[2]
drugs_kept = names(which(na_rows/num_cells <0.5))

#quantify reduction
length(drugs_kept)/dim(cancer.data)[1]

#reduce
cancer.data = cancer.data[drugs_kept, ]
sum(is.na(cancer.data))/(dim(cancer.data)[1]*dim(cancer.data)[2])

library(zoo)
cancer.data = t(cancer.data)
cancer.data = na.aggregate(cancer.data)
cancer.data = t(cancer.data)

### PCA
cancer.cov = cov(cancer.data)
cancer.eigen = eigen(cancer.cov)

#loadings hard to make sense of (could average over cancer type and plot something smaller...)
x11()
par(mfrow=c(5,1))
for(i in 1:5){
  barplot(cancer.eigen$vectors[,i])
}

#visualization
x11()
par(mfrow=c(1,2))
plot(cancer.eigen$values, xlab="pc", ylab="variance")
plot(cumsum(cancer.eigen$values)/sum(cancer.eigen$values), xlab="pc", ylab="explained variance", type='b')
abline(h=0.85, col='red', lty=2)


thresh = 0.9
n.comp = min(which(cumsum(cancer.eigen$values)/sum(cancer.eigen$values)>=thresh))

dev.off()

cancer.scores = cancer.data%*%cancer.eigen$vectors[,1:n.comp]
dim(cancer.scores)

#now let's cluster and generate the curve (use kmeans for SSres and such )
J = c(NULL)

for(i in 2:20){
  cluster.obj = kmeans(cancer.scores, i)
  J = c(J,cluster.obj$tot.withinss/(cluster.obj$tot.withinss + cluster.obj$betweenss))
}

plot(J, type = 'b')
abline(v=6, col='red', lty=2)

n.clusts = 6
cluster.obj = kmeans(cancer.scores, n.clusts)

#let's guess like 8
#visualize 
x11()
plot(cancer.scores[,1:2], pch = 19, col = cluster.obj$cluster, xlab = "PC1", ylab = "PC2")

#we should check the gaussian properties of these bad boys 
load("/Users/nathanielnethercott/Desktop/School/Polimi/2021:2022/AS/Labs-20220417/Lab 5/mcshapiro.test.RData")

#NEED TO VERIFY GAUSSIANITY
mcshapiro.test(cancer.scores[which(cluster.obj$cluster == 1), ])

#let's compare with our actual labels now 
drug.labels = read.delim("../drugs/drugs_labelled.csv", sep = ",", row.names = 1)
dim(drug.labels)
drug.labels = drug.labels[which(drug.labels$NAME%in%rownames(cancer.data)),]

#fuck i need to resolve the name 1/2 vs name-2 convention 


drug.clusters = cluster.obj$cluster

#main table for counts of labelled drugs
ref = table(drug.labels$CancerType)

for(i in 1:n.clusts){
  clust_drugs = names(which(drug.clusters == i))
  clust_cancer_table = table(drug.labels[which(drug.labels$NAME%in%clust_drugs),"CancerType"])
  
  print(paste("CLUSTER no. ", i))
  print(clust_cancer_table)
}



#IDEAS:
#' - higher representation of one cancer might result in us omitting one or low samples from
#'   other cancers which are not as prevalent when we compute the PCA (it is a linear combo of cancers)
#' 

#ISSUES:
#' - we can't easily test gaussianity for the clusters using mcshapiro test util since p>n 
#' 


#----------------------------------------------------------------------
#try the same process again but with new `features` for each drug 
rm(list=ls())

source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")

#borrowing syntax from `auc_exp4.R`
cancers = cancer_types
cancer.data = block_dat(cancers, auc)

na_rows = sort(rowSums(is.na(cancer.data)), decreasing=TRUE)
num_cells = dim(cancer.data)[2]
drugs_kept = names(which(na_rows/num_cells <0.5))

#reduce
cancer.data = cancer.data[drugs_kept, ]

library(zoo)
cancer.data = t(cancer.data)
cancer.data = na.aggregate(cancer.data)
cancer.data = t(cancer.data)

#this is gonna be bricked 
reduced = data.frame(rowMeans(block_dat(c(cancers[1]), cancer.data)))
for(i in 2:length(cancers)){
  reduced = cbind(reduced, rowMeans(block_dat(c(cancers[i]), cancer.data)))
}

colnames(reduced) = cancers
View(reduced)

#we decrease info to 4% of original 
dim(reduced)[2]/dim(cancer.data)[2]

############ PCA ########################
reduced.pc = princomp(reduced, scores=TRUE)

x11()
par(mfrow=c(1,2))
plot(reduced.pc$sdev, xlab="pc", ylab="variance")
plot(cumsum(reduced.pc$sdev)/sum(reduced.pc$sdev), xlab="pc", ylab="explained variance", type='b')
abline(h=0.85, col='red', lty=2)

thresh = 0.9
n.comp = min(which(cumsum(reduced.pc$sdev)/sum(reduced.pc$sdev)>=thresh))

#try the clustering w cost function again 
reduced.scores = reduced.pc$scores[,1:n.comp]
dim(reduced.scores)

#now let's cluster and generate the curve (use kmeans for SSres and such )
J = c(NULL)

for(i in 2:20){
  cluster.obj = kmeans(reduced.scores, i)
  J = c(J,cluster.obj$tot.withinss/(cluster.obj$tot.withinss + cluster.obj$betweenss))
}

plot(J, type = 'b')
abline(v=6, col='red', lty=2)

n.clusts = 6
cluster.obj = kmeans(reduced.scores, n.clusts)

#plot 
plot(reduced.scores[,1], reduced.scores[,2], col = cluster.obj$cluster, pch = 19)

#--------------------------------------------
#let's compare with our actual labels now 
drug.labels = read.delim("../drugs/drugs_labelled.csv", sep = ",", row.names = 1)
dim(drug.labels)
drug.labels = drug.labels[which(drug.labels$NAME%in%rownames(reduced)),]

drug.clusters = cluster.obj$cluster

#main table for counts of labelled drugs
ref = table(drug.labels$CancerType)

for(i in 1:n.clusts){
  clust_drugs = names(which(drug.clusters == i))
  clust_cancer_table = table(drug.labels[which(drug.labels$NAME%in%clust_drugs),"CancerType"])
  
  print(paste("CLUSTER no. ", i))
  print(clust_cancer_table)
}
