######################## IDEAS DISCUSSED ########################
# - group effects introduced by collection lab on auc
# - 

#################################################################

features = read.delim(file.path(getwd(), "physiological.csv"), header = TRUE, sep=",")
names(features)

#let's look at a boxplot of auc variance based on the `CELL_LINE-SOURCE`
#manually add this instead of modifying the original constructed dataset
  #- might only make sense to consider specific cancers with this in mind (e.g. all breast but different sources)

path_to_data = "/Users/nathanielnethercott/Desktop/School/Polimi/2021:2022/AS/AS_Project_2022/Dataset/"
samples = read.delim(file.path(path_to_data, "data_clinical_sample.txt"), header=TRUE, comment.char = '#')

#really hoping this preserves the order 
head(samples$PATIENT_ID[samples$SAMPLE_ID%in%features$PATIENT_ID])
head(features$PATIENT_ID)

features$CELL_LINE_SOURCE = samples$CELL_LINE_SOURCE[samples$SAMPLE_ID%in%features$PATIENT_ID]

x11()
barplot(table(features$CELL_LINE_SOURCE))

#get rid of empty sex values
features[features == ""] = NA

sort(table(features$CELL_LINE_SOURCE), decreasing = TRUE)
sort(table(features$CANCER_TYPE), decreasing = TRUE)


#----------------------------------------
#consider a subset of cells and consider the same type of plot 
source("../utils/nate_utils.R")

#select cancer cells from main dataset
cancer_name = "Breast Cancer"
cancer.dat = features[features$CANCER_TYPE == cancer_name, ]

#----------------------------------------
#' ADDITIONAL STEP IN PIPELINE FOR RECOMPUTING AUC 
#' some logic here where we redefine the auc to be an average over relevant drugs 
#' bleh 

drug.labels = read.delim("../drugs/drugs_labelled.csv", sep = ",", row.names = 1)
dim(drug.labels)

#okay this is where you can use the string utility I made (i was stupid haha)
temp = features[features$CANCER_TYPE == cancer_name, ]
cancer_type_nate = unique(str_util(temp$PATIENT_ID))
rm(temp)

#now we can cross-reference the list of known drugs for treating that cancer 
drugs.relevant = drug.labels$NAME[which(drug.labels$CancerType == cancer_type_nate)]

#handle the auc import 
auc = read.delim(file.path(path_to_data, "data_drug_treatment_auc.txt"), header=TRUE, comment.char = '#')
auc.new = auc[,-c(1,3,4)]
rm(auc)
auc.new = auc.new[which(auc.new$NAME%in%drugs.relevant),]
rownames(auc.new) = auc.new$NAME
auc.new = auc.new[,which(colnames(auc.new)%in%cancer.dat$PATIENT_ID)]

#recompute the auc across these new drugs now 
auc.mean = colMeans(auc.new, na.rm=TRUE)

#reassign the auc column
for(cell in cancer.dat$PATIENT_ID){
  cancer.dat[which(cancer.dat$PATIENT_ID == cell),"AUC"] = auc.mean[cell]
}

#remove any fully NA samples
cancer.dat = cancer.dat[!is.na(cancer.dat$AUC),]

#----------------------------------------

#fit some anova models to quantify relationships?

#BEFORE DOING THIS WE NEED 
  #- guassianity w.r.t each group 
  #- homoschedasicity 

#need to omit some rows pertaining to missing categorical variable (e.g. in sex case, but not ethnicity)

#----------------------------------------
#AN EXTENSIVE STUDY OF THIS CANCER WE CHOSE 
#age vs auc (for predictive reasons)
plot(cancer.dat$AGE, cancer.dat$AUC, main = cancer_name, xlab = 'age', ylab='auc', pch = 19)

#features pairwise, look for apparent correlations 
names(cancer.dat)
x11()
pairs(cancer.dat[,c(3,5,6,9)], pch = 19)
dev.off()

#categorical 
table(cancer.dat$SEX)
table(cancer.dat$ETHNICITY)
table(cancer.dat$CELL_LINE_SOURCE)


#1. CELL_LINE_SOURCE
#considering if the cell line source has an effect on the regenerative properties of a cell 
temp.dat = cancer.dat
sum(is.na(temp.dat))

#qualitative w less than 3 samples (so we can use R's shapiro test)
#cell line source 
labs_to_keep = names(which(table(temp.dat$CELL_LINE_SOURCE)>3))
labs_to_keep
temp.dat = temp.dat[which(temp.dat$CELL_LINE_SOURCE%in%labs_to_keep),]

#check our pre-reqs for applying an aov analysis 
attach(temp.dat)

#normality 
tapply(AUC, CELL_LINE_SOURCE, shapiro.test)
#homoschedasticity 
bartlett.test(AUC, CELL_LINE_SOURCE)

detach(temp.dat)

fit1 = aov(AUC~CELL_LINE_SOURCE, data = temp.dat)
summary.aov(fit1)

#visualize
x11()
boxplot(cancer.dat$AUC~cancer.dat$CELL_LINE_SOURCE, xlab='lab', ylab='auc', col='gold')

#[NATE] manually rederiving the f-test for aov group significance 
group.means = tapply(temp.dat$AUC, temp.dat$CELL_LINE_SOURCE, mean)
pop.mean = mean(temp.dat$AUC)
g = length(unique(temp.dat$CELL_LINE_SOURCE))
n = dim(temp.dat)[1]

group.counts = table(temp.dat$CELL_LINE_SOURCE)
mean.var = (1/(g-1))*sum(group.counts*(group.means - rep(pop.mean, g))^2)
group.var = (1/(n-g))*sum(fit1$residuals^2)

1 - pf(mean.var/group.var, g-1, n-g)
summary.aov(fit1)
#[NATE] okay that worked 


#2. ETHNICITY 
temp.dat = cancer.dat #reassign

#ethnicity 
eth_to_keep = names(which(table(temp.dat$ETHNICITY)>3))
eth_to_keep
temp.dat = temp.dat[which(temp.dat$ETHNICITY%in%eth_to_keep),]

#check our pre-reqs for applying an aov analysis 
attach(temp.dat)

#normality 
tapply(AUC, ETHNICITY, shapiro.test)
#homoschedasticity 
bartlett.test(AUC, ETHNICITY)

detach(temp.dat)

boxplot(temp.dat$AUC~temp.dat$ETHNICITY, xlab = 'ethnicity', ylab = 'auc', col='gold')
#doesn't appear like there's gonna be a convincing result 

fit2 = aov(AUC~ETHNICITY, data=temp.dat)
summary.aov(fit2)


#3. SEX 
temp.dat = cancer.dat #reassign

sex_to_keep = names(which(table(temp.dat$SEX)>3))
sex_to_keep
temp.dat = temp.dat[which(temp.dat$SEX%in%sex_to_keep),]

attach(temp.dat)

#normality 
tapply(AUC, SEX, shapiro.test)
#homoschedasticity 
bartlett.test(AUC, SEX)

detach(temp.dat)

boxplot(temp.dat$AUC~temp.dat$SEX, xlab = 'sex', ylab = 'auc', col='gold')

x11()
par(mfrow=c(1,2))
hist(temp.dat[which(temp.dat$SEX == "Male"),]$AUC, main = "Male", xlab = 'auc')
hist(temp.dat[which(temp.dat$SEX != "Male"),]$AUC, main = "Female", xlab = 'auc')
dev.off()

fit3 = aov(AUC~SEX, data=temp.dat)
summary.aov(fit3)




