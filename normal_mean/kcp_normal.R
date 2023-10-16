require(gtools)
library(boot)
library(magrittr)
library(MASS)
library(Hotelling)
library(gTests)
library(matrixStats)
library(ade4)
library(mvtnorm)

##### new package
require("gSeg")
require(Rcpp)
require("kerSeg") 
source("objTest_fctns.R")
source("depth_CPD_func.R") 
source("ecp_distmat_input.R")
source("kcp_distmat_input.R")
sourceCpp('energyChangePoint.cpp')

delta<-as.numeric(commandArgs(TRUE)[1])
monte_carlo<-as.numeric(commandArgs(TRUE)[2])
# for mean_diff case
S1<-load("Sigma_30.Rdata")
S1<-get(S1)
S2<-load("Sigma_90.Rdata")
S2<-get(S2)
S3<-load("Sigma_180.Rdata")
S3<-get(S3)

n1=100;n2=200;n=n1+n2
num_permut<-1000
result<-list()
for (i in seq_along(1:3)){
  
  # There are four different cases and we have different parameters for them.
  # so we need to change the corresponding data generation method in the code and also input parameters, i.e., this .R file and parameter.sh file 
  
  # 1. mean_diff: delta in seq(0,1,0.1), 
  # 2. scale_diff: delta in seq(0,0.4,0.04), 
  # 3. mixture gaussian: delta in seq(0,1,0.1), 
  # 4. heavy tail: delta in seq(2,22,2)
  
  
  p=c(30,90,180)[i] # for mean_diff, scale_diff, mixture gaussian cases
  #p=c(5,15,60)[i] # for heavy tail case
  set.seed(10000*monte_carlo+100*delta)
  I<-diag(x = 1, p, p)
  #### 1. mean diff
  Sigma=list(S1,S2,S3)[[i]]
  Data<-rbind(mvrnorm(n1,mu=rep(0,p),Sigma=Sigma),mvrnorm(n2,mu=c(rep(delta,p)),Sigma=Sigma))
  
  #### 2. scale diff
  # Data<-rbind(mvrnorm(n1,mu=rep(0,p),Sigma=0.8*I),mvrnorm(n2,mu=c(rep(0,p)),Sigma=(0.8-delta)*I))
  
  
  #### 3. mixture gaussian distributions
  #A<-rbinom(n2,1,0.5)
  #mu<-c(rep(delta,0.1*p),rep(0,0.9*p))
  #Z1<-mvrnorm(n2,-mu,I);Z2<-mvrnorm(n2,mu,I)
  #Data<-rbind(mvrnorm(n1,mu=rep(0,p),Sigma=I),A*Z1+(1-A)*Z2)
  
  
  #### 4. heavy tailed distribution: different degrees of freedom
  # Data<-rbind(mvrnorm(n1,mu=rep(0,p),Sigma=I),matrix(rt(n2*p, df=delta),nrow=n2))
  
  distmat<-as.matrix(dist( Data, method = 'euclidean' ))
  depth_result<-depth_CPD(distmat,num_permut =1,c=0.1)
  
  #### Graph CPD
  E1 = mstree(dist( Data, method = 'euclidean' ),ngmax = 5)
  result_graph = gseg1(nrow(distmat),E1, statistics="all")
  #### Energy CPD
  result_ecp<-e.divisive_distmat(D=distmat,sig.lvl=.05,R=999,k=NULL,min.size=30,alpha=1)
  
  #### kcp
  result_kcp<-kcpa_distmat(D=distmat,L=nrow(distmat)/sqrt(log(nrow(distmat))))
  ####kerSeg
  K = gaussiankernel(Data) # Gaussian kernel matrix
  result_kerseg = kerseg1(nrow(distmat), K, pval.perm=TRUE, B=1000)
  
  r<-list(depth_result,result_graph,result_ecp,result_kcp,result_kerseg)
  names(r)<-c("depth",'graph','ecp','kcp','kerseg')
  result[[i]]<-r
}

path<-paste("delta_",delta,'_run_',monte_carlo,'.Rdata',sep="")
save(result, file=path)


