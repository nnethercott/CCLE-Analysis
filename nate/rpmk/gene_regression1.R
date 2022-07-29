######################## IDEAS DISCUSSED ########################
#' - first find relevant genes which characterize breast cancer through regression 
#'  tree best matching real label accuracy 
#' - look at google scholar results to find drugs relevant to breast cancer 
#' **in the worst case we just use some of the genes from online or whatever 
#################################################################

source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")

#using my rpkm dataset, not the ~exact~ same one as luca 
View(rpkm)

#################################################################
######################### SETUP DATA #############################
#################################################################
M <- scale(rpkm, center = T, scale = T) 
threshold_l <- 1
threshold_h <- 100

row_var = apply(M, 1, var) #apply over rows: variability along the genes

x11()
par(mfrow=c(1,2))
#plot 1
plot((row_var))
abline(h=threshold_l, col='red')
abline(h=threshold_h, col='red')
#plot 2
plot(row_var[row_var > threshold_l & row_var < threshold_h])

#keep only rows above the lower threshold for original rpkm
rpkm.clean = rpkm[row_var > threshold_l,]


#~~~~~~~~~~~~~~~~~~~~ (extending stuff) ~~~~~~~~~~~~~~~~~~~~~~~~~~
cancers = c("BREAST")
#probably add a step here where we create all tuples from cancers of interest 

cancers.rpkm = block_dat(cancers, rpkm.clean)
#View(cancers.rpkm)
dim(cancers.rpkm)

cancers.auc = block_dat(cancers, auc)
dim(cancers.auc)

# Target variable
# y <-as.factor(labels2.real)
#let's create a tier system for the severity of cancers across relevant drugs 
drug.labels = read.delim("../drugs/drugs_labelled.csv", sep = ",", row.names = 1)
drug.labels = drug.labels[which(drug.labels$CancerType!=""),]

drugs_relevant = c(NULL)
for(i in 1:dim(drug.labels)[1]){
  for(cancer in cancers){
    if(cancer%in%unlist(strsplit(drug.labels$CancerType[i], ", "))){
      drugs_relevant = c(drugs_relevant, drug.labels$NAME[i])
    }
  }
  
}
drugs_relevant

#now get a relevant auc average 
auc_relevant = cancers.auc[rownames(cancers.auc)%in%drugs_relevant,]
auc_relevant = t(auc_relevant)

#get rid of na's 
library(zoo)
auc_relevant = na.aggregate(auc_relevant)
y_ = rowMeans(auc_relevant)

#transform this bad boy maybe... ***************

#okay now create the tiers, first need midpoints to which we minimize dist 
n_tiers = 2
step = (max(y_)-min(y_))/(n_tiers)
midpoints_ = seq(from=min(y_),to=max(y_),by=step)
midpoints_

midpoints = c((midpoints_[1]+midpoints_[2])/2)
for(i in 2:(length(midpoints_)-1)){
  midpoints = c(midpoints, (midpoints_[i]+midpoints_[i+1])/2)
}
midpoints

#bin the data 
library(purrr)
y = map(y_, function(x) which.min(abs(x-midpoints)))
y = as.numeric(y)

#visualize
plot(y_,y, pch=19)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################
######################### FEATURE SELECTION #############################
library("dplyr")
data_expression <- t(cancers.rpkm)
data_expression <- data.frame(data_expression)

x <- as.matrix(data_expression)
response = y_

############################# lasso regression ##########################
library(glmnet) #lasso

fit.lasso = cv.glmnet(x,response, alpha=1)
best_lambda <- fit.lasso$lambda.min
best_lambda

plot(fit.lasso) 

best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
names(best_model)
lasso.coef = coef(best_model)

sort(abs(lasso.coef), decreasing=TRUE)

threshold = 1e-7
indices = which(abs(lasso.coef)>threshold)
indices
genes.lasso = colnames(data_expression)[indices[-1]-1]
genes.lasso

#explore some stuff w lasso, cumulative effect of the coefficients 
plot(seq(1:length(lasso.coef)), cumsum(sort(abs(lasso.coef), decreasing=TRUE))/sum(abs(lasso.coef)))

#what are the top few biggest ones 
order(abs(lasso.coef), decreasing = TRUE)[1:6]
genes.lasso_reduced = colnames(data_expression)[c(305, 425, 260,  59, 409)-1]

#################################################################
######################### PREDICTIVE ############################
#################################################################

#first see if the original lasso regression model makes any sense 
auc_predicted <- predict(best_model, s = best_lambda, newx = x)
auc_predicted[1:10] #this is ass (since larger than 1)

target_genome = data_expression[,which(colnames(data_expression)%in%genes.lasso)]
target_genome_reduced = data_expression[,which(colnames(data_expression)%in%genes.lasso_reduced)]

#train a logistic regression model on our new subset of genes??
#SUPER BIG ISSUE IS THAT GENE EXPRESSION IS LIKE ZERO FOR SOME PEOPLE ON `IMPORTANT` GENES
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(target_genome_reduced[,i], response, pch=19, xlab = genes.lasso_reduced[i])
}

x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(log(target_genome_reduced[,i]), response, pch=19, xlab = paste("log(", genes.lasso_reduced[i],")"))
}


#this should be the logistic regression setup 
x11()
par(mfrow=c(2,2))
plot(log(target_genome_reduced[,1]), y-1, pch=19)
plot(log(target_genome_reduced[,2]), y-1, pch=19)
plot(log(target_genome_reduced[,3]), y-1, pch=19)
plot(log(target_genome_reduced[,4]), y-1, pch=19)

training = data.frame(g1 = target_genome_reduced[,1], 
                      g2 = target_genome_reduced[,2],
                      g3 = target_genome_reduced[,3],
                      g4 = target_genome_reduced[,4], 
                      response = y-1)


logistic.model = glm(response~
                       g1+g2+g3+g4,
                     data = training, family=binomial(link=logit))


#################### INCORPORATE PHYSIOLOGICAL #############################
#as was done before... 
features = read.delim(file.path(getwd(), "../physiological/physiological.csv"), header = TRUE, sep=",")
names(features)
path_to_data = "/Users/nathanielnethercott/Desktop/School/Polimi/2021:2022/AS/AS_Project_2022/Dataset/"
samples = read.delim(file.path(path_to_data, "data_clinical_sample.txt"), header=TRUE, comment.char = '#')
features$CELL_LINE_SOURCE = samples$CELL_LINE_SOURCE[samples$SAMPLE_ID%in%features$PATIENT_ID]
features$AGE_ENCODED = as.numeric(factor(cut(features$AGE, breaks = c(0,25,50,75,100))))-1

cancer.features = features[which(features$PATIENT_ID%in%rownames(target_genome)),]



####### SOME PLOTS WITH COLOR MAPS AS A FINAL PRAYER FOR RESULTS ############
table(cancer.features$ETHNICITY)
table(cancer.features$SEX)
table(cancer.features$CELL_LINE_SOURCE)
table(cancer.features$AGE_ENCODED)

#ethnicity 
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(log(target_genome_reduced[,i]), response, pch=19, xlab = paste("log(", genes.lasso_reduced[i],")"),
       col = as.numeric(factor(cancer.features$ETHNICITY))+1)
}

#sex 
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(log(target_genome_reduced[,i]), response, pch=19, xlab = paste("log(", genes.lasso_reduced[i],")"),
       col = as.numeric(factor(cancer.features$SEX))+1)
}


#cell line 
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(log(target_genome_reduced[,i]), response, pch=19, xlab = paste("log(", genes.lasso_reduced[i],")"),
       col = as.numeric(factor(cancer.features$CELL_LINE_SOURCE))+1)
}

#age 
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(log(target_genome_reduced[,i]), response, pch=19, xlab = paste("log(", genes.lasso_reduced[i],")"),
       col = as.numeric(factor(cancer.features$AGE_ENCODED))+1)
}



############## not logged ############
#ethnicity 
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(target_genome_reduced[,i], response, pch=19, xlab = genes.lasso_reduced[i],
       col = as.numeric(factor(cancer.features$ETHNICITY))+1)
}

#sex 
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(target_genome_reduced[,i], response, pch=19, xlab = genes.lasso_reduced[i],
       col = as.numeric(factor(cancer.features$SEX))+1)
}


#cell line 
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(target_genome_reduced[,i], response, pch=19, xlab = genes.lasso_reduced[i],
       col = as.numeric(factor(cancer.features$CELL_LINE_SOURCE))+1)
}


#age 
x11()
par(mfrow=c(2,2))
for(i in 1:4){
  plot(target_genome_reduced[,i], response, pch=19, xlab =  genes.lasso_reduced[i],
       col = as.numeric(factor(cancer.features$AGE_ENCODED))+1)
}


########## 3d stuff ########
library(rgl)
par3d(windowRect = c(20, 30, 800, 800))
plot3d(target_genome_reduced[,1], target_genome_reduced[,2], response, pch=19, col='blue', size=10,
       xlab = paste("log(",genes.lasso_reduced[1],")"), ylab=paste("log(",genes.lasso_reduced[2],")"))

####################################################################################
############################### LINEAR REGRESSION ##################################
#full model (much better)
W = data.frame(cbind(target_genome,y=response))
fit1 = lm(y~., data=W)
summary(fit1)

names(fit1)
sort(abs(coefficients(fit1)),decreasing=T)

qqnorm(fit1$residuals)
qqline(fit1$residuals)

#training dataset 
training = W[1:35,]
test = W[36:46,]

fit2 = lm(y~., data=training)
summary(fit2)

auc_pred = predict(fit2, test)
sqrt(sum((test$y-auc_pred)^2)/length(test$y))

#visualize
plot(seq(1:length(test$y)), auc_pred, col = 'blue', pch=19, ylab="auc", xlab="sample")
points(seq(1:length(test$y)), test$y, col='red',pch=19)


#try another method considering age and stuff as dummy variables 
W2 = data.frame(cbind(genes = target_genome,eth=as.numeric(factor(cancer.features$ETHNICITY))-1, y=response))
fit2 = lm(y~., data=W2)
summary(fit2)


############ REDUCED #############
Z = data.frame(cbind(target_genome_reduced, y=response))
fit3 = lm(y~., data=Z)
summary(fit3)

Z1 = data.frame(cbind(genes=target_genome_reduced,age=as.numeric(factor(cancer.features$CELL_LINE_SOURCE))-1, y=response))
fit4 = lm(y~age + genes.CCT7+genes.LDHA+genes.AC107016.2+genes.RPL23A+genes.CD79B+
            genes.CCT7:age+genes.LDHA:age+genes.AC107016.2:age+genes.RPL23A:age+genes.CD79B:age, 
            data=Z1)
summary(fit4)

#####################################################################################
#LUCA THING THAT SOUNDS FAKE ########################################################

library(ISLR)
library(leaps)

W = data.frame(cbind(target_genome,y=response))

p=13
regfit.full <- regsubsets(y~., data=W, nvmax=p)
summary(regfit.full)

reg.summary <- summary(regfit.full)

x11()
par(mfrow=c(1,3))
plot(reg.summary$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
abline(v=which.max(reg.summary$adjr2), col = "red", lty =2)
plot(reg.summary$rss,xlab="Number of Variables",ylab="RSS",type="b")

#extract relevant genes
genes.regfit = names(coef(regfit.full,10))[-1]
regfit_target_genome = target_genome[,colnames(target_genome)%in%genes.regfit]

W_regfit= data.frame(cbind(regfit_target_genome,y=response))
fit.exhaustive = lm(y~., data=W_regfit)
summary(fit.exhaustive)

qqnorm(fit.exhaustive$residuals)
qqline(fit.exhaustive$residuals)

#try predicting 
plot(seq(1:length(fit.exhaustive$fitted.values)),fit.exhaustive$fitted.values, pch=19, col='red',ylim=c(0.6,1))
points(seq(1:length(fit.exhaustive$fitted.values)), auc_relevant[,6],col='black')


#now do the combined version #############
#load luca genes
load("./selected_genes_using_breast_nervous.Rdata")
genes.luca = good_genes[1:13]
genes.luca = c(genes.luca, c("CSTL", "RP11-67L3.5", "DHCR24"))

genes.combined = c(genes.regfit, genes.luca)
genome = data_expression[,which(colnames(data_expression)%in%genes.combined)]

W.combined = data.frame(cbind(genome,y=response))

p=13
regfit.combined <- regsubsets(y~., data=W.combined, nvmax=p)
summary(regfit.combined)

reg.summary_combined <- summary(regfit.full)

x11()
par(mfrow=c(1,3))
plot(reg.summary_combined$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(reg.summary_combined$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
abline(v=which.max(reg.summary_combined$adjr2), col = "red", lty =2)
plot(reg.summary_combined$rss,xlab="Number of Variables",ylab="RSS",type="b")

#extract relevant genes
genes.regfit_combined = names(coef(regfit.combined,8))[-1]
regfit_target_genome_combined = genome[,colnames(genome)%in%genes.regfit_combined]

W_regfit_combined= data.frame(cbind(regfit_target_genome_combined,y=response))
fit.exhaustive_combined = lm(y~., data=W_regfit_combined)
summary(fit.exhaustive_combined)

x11()

qqnorm(fit.exhaustive_combined$residualssa)
qqline(fit.exhaustive_combined$residuals)

#see something 
genes.luca%in%genes.regfit_combined


#consider cross-validation on the linear model regression the chosen genes
library(caret)
k_ = 5
flds <- createFolds(response, k = k_, list = TRUE, returnTrain = FALSE)

cv_rsq = c(NULL)
cv_mse = c(NULL)

for(i in 1:k_){
  fold = as.numeric(unlist(flds[i]))
  train = W_regfit_combined[-fold,]
  test = W_regfit_combined[fold,]
  
  fit.temp = lm(y~., data=train)
  cv_rsq = c(cv_rsq, summary(fit.temp)$adj.r.squared)
  
  #predicted 
  pred = predict(fit.temp, test)
  cv_mse = c(cv_mse, sqrt(sum((test$y - pred)^2)/length(test$y)))
}
cv_rsq
sqrt(cv_mse)
