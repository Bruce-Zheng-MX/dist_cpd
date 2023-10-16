library(igraph)
require("distrEx")
library(mvtnorm)
require(gtools)
library(boot)
library(magrittr)
library(MASS)
library(Hotelling)
library(gTests)
library(matrixStats)
library(ade4)


MMD_test<-function(D){
  
  bw<-median(D)
  gram_matrix<-exp(D^2*(-1/(2*bw^2))) 
  T<-nrow(D)
  
  obj_vals<-matrix(0,T)
  for (t in seq(2,T-1)){
    term1 = (T-t)/(t*T)*sum(gram_matrix[1:t, 1:t])
    term2 = -2/T*sum(gram_matrix[1:t, t:T])
    term3 = t/((T-t)*T)*sum(gram_matrix[t:T, t:T])
    obj_vals[t] = term1 + term2 + term3
  }
  return(obj_vals)
}


##### new package
require("gSeg")
require(Rcpp)
require("kerSeg") 
source("objTest_fctns.R")
source("depth_CPD_func.R") 
source("ecp_distmat_input.R")
source("kcp_distmat_input.R")
sourceCpp('energyChangePoint.cpp')

n1=100;n2=200;n=n1+n2
num_permut<-1000
v1<-0.0000000001
delta<-as.numeric(commandArgs(TRUE)[1])
monte_carlo<-as.numeric(commandArgs(TRUE)[2])

Data<-c()
for (i in 1:n1){
  g<-as.matrix(laplacian_matrix(sample_pa(200,power=v1,directed=FALSE)))
  Data<-rbind(Data,c(2*g[upper.tri(g)],diag(g)))
}
for (i in 1:n2){
  g<-as.matrix(laplacian_matrix(sample_pa(200,power=delta,directed=FALSE)))
  Data<-rbind(Data,c(2*g[upper.tri(g)],diag(g)))
}

distmat<-as.matrix(dist( Data, method = 'manhattan' ))
depth_result<-depth_CPD(distmat,num_permut =num_permut,c=0.1)

#### Graph CPD
E1 = mstree(dist( Data, method = 'manhattan' ),ngmax = 5)
result_graph = gseg1(nrow(distmat),E1, statistics="all")
#### Energy CPD
result_ecp<-e.divisive_distmat(D=distmat,sig.lvl=.05,R=999,k=NULL,min.size=30,alpha=1)



#### MMD test
MMD<-c()
for (j in 0:num_permut){
  if (j!=0){
    set.seed(j)
    ind<-sample(nrow(distmat))
    D<-distmat[ind,ind]
  }else{
    D<-distmat
  }
 result_mmd<-MMD_test(D)
 loc<-which.max(result_mmd)
 ob_stat<-result_mmd[loc]
 MMD<-rbind(MMD,c(loc,ob_stat))
}
  
r<-list(depth_result,result_graph,result_ecp,MMD)
names(r)<-c("depth",'graph','ecp','mmd')



path<-paste("delta_",delta,'_run_',monte_carlo,'.Rdata',sep="")
save(r, file=path)

