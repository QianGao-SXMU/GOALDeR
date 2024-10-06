# GOALDeR: Generalized Outcome-Adaptive LASSO and Doubly Robust Estimate  
Implementation of the Generalized Outcome-Adaptive LASSO and Doubly Robust Estimate (GOALDeR) and Doubly robust estimator using super learner to fit the GPS model and the outcome model (SL-DR).  
# Installation  
library(devtools)  
install_local(“/path/to/lqa_1.0-3.tar.gz”)  
install_github("heike/extracat")  
install_github("QianGao-SXMU/GOALDeR")  
# Example  
## Simulation Data  
n<-500  
p<-200  
mean_x <- 0  
sig_x <- 1  
rho <- 0.2  
Beta<-0.7  
var.list<-paste("X",1:p,sep="")  
Sigma <- matrix(rho\*sig_x^2,nrow=length(var.list),ncol=length(var.list))  
diag(Sigma) <- sig_x^2  
Mean_x <- rep(mean_x,p)  
set.seed(123)  
Data <- as.data.frame(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma,empirical = FALSE))  
colnames(Data)<-var.list  
betaOut<-c(1,1,1,1,0,0)  
betaTrt<-c(0.5,0.5,0,0,1,1)  
TrueVar<-as.matrix(Data[,c(1:6)])  
Data$Trt<-TrueVar%\*%betaTrt+rnorm(n=n,mean=0,sd=1)  
Data$Y<-Beta\*Data$Trt+TrueVar%\*%betaOut+rnorm(n=n,mean=0,sd=1)  
# GOALDeR and SL-DR  
library(Matrix)  
library(MASS)  
library(lqa)  
library(extracat)  
library(survey)  
library(cdcsis)  
library(extracat)  
library(independenceWeights)  
GOALDeR(data=Data,covar=NULL,var.list=var.list,Trt="Trt",out="Y",gamma_convergence=2)  
SL_DR(Y=Data$Y,A=Data$Trt,X=Data[,var.list],treat_mod = "Normal")
