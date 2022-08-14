library(BGLR)
library(plyr)
library(data.table)

#Calculating program execution time
Beginning_time<-'t1'
Ending_time<-'t2'
t1<-proc.time()

#Loading dataset
#Loading wheat dataset
load("D:/Bayesian models/standarizedData_univariate.RData")
wheat_pheno<-pheno
wheat_gene<-X
#Importing the German Holstein Dairy cattle Dataset
cattle_pheno<-read.table(file="FileS2.txt",header = TRUE)#phenotypic data
cattle_gene<-read.table(file="cattle_genotypes.txt",header = TRUE)#genomic data

#parameters
seed=123#set the seed for the random number generator
nIter<-50000 
burnIn<-25000

#10-fold cross validation
folds<-10
y<-wheat_pheno[,c(6)]#Read trait (GP) in the sixth column
set.seed(seed) #Set seed for the random number generator
sets<-rep(1:10,200)[-1]
sets<-sets[order(runif(nrow(wheat_gene)))]
COR.CV<-rep(NA,times=(folds+1))
MSE.CV<-rep(NA,times=(folds+1))
MAE.CV<-rep(NA,times=(folds+1))
RMSE.CV<-rep(NA,times=(folds+1))
w<-rep(1/nrow(wheat_gene),folds) ## weights for pooled correlations,MAE,RMSE,and MSE
yHatCV<-numeric()
for(fold in 1:folds)
{
  yNa<-y
  whichNa<-which(sets==fold)
  yNa[whichNa]<-NA
  ETA<-list(list(X=wheat_gene,model='BRR'))#other models:'BL','BayesA','BayesB','BayesC'
  fm<-BGLR(y=yNa,ETA=ETA,
           nIter=nIter,burnIn=burnIn)
  yHatCV[whichNa]<-fm$yHat[fm$whichNa]
  w[fold]<-w[fold]*length(fm$whichNa)
  COR.CV[fold]<-cor(fm$yHat[fm$whichNa],y[whichNa])
  MAE.CV[fold]= mean(abs(fm$yHat[fm$whichNa]-y[whichNa]))
  MSE.CV[fold]<-mean((fm$yHat[fm$whichNa]-y[whichNa])^2)
  RMSE.CV[fold]<-sqrt(mean((fm$yHat[fm$whichNa]-y[whichNa])^2))
  real_value<-y[whichNa]
  pred_value<-fm$yHat[fm$whichNa]
  write.csv(real_value,"BRRrealPHE6_CV10.csv")
  write.csv(pred_value,"BRRpredPHE6_CV10.csv")
}

COR.CV[11]<-mean(COR.CV[1:10])
MAE.CV[11]<-mean(MAE.CV[1:10])
MSE.CV[11]<-mean(MSE.CV[1:10])
RMSE.CV[11]<-mean(RMSE.CV[1:10])
COR.CV
MAE.CV
MSE.CV
RMSE.CV

#Calculating program execution time
t2 = proc.time()
running_time = t2 - t1
print(paste0('running time£º',running_time[3][[1]],'s'))



