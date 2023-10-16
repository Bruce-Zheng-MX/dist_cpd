require(gtools)
library(boot)
library(magrittr)
library(MASS)
#library(energy)
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

delta<-as.numeric(commandArgs(TRUE)[1])
monte_carlo<-as.numeric(commandArgs(TRUE)[2])

n1=100;n2=200;n=n1+n2
num_permut<-1000
I<-diag(1,2)

# mean shift: delta in seq(0,1,0.1), 
z1<-mvrnorm(n=n1,mu=c(0,0),Sigma = 0.25*I)
z2<-mvrnorm(n=n2,mu=c(delta,0),Sigma = 0.25*I)

#scale shift: delta in seq(0,0.4,0.04), 
#z1<-mvrnorm(n=n1,mu=c(0,0),Sigma = 0.16*I)
#z2<-mvrnorm(n=n2,mu=c(0,0),Sigma = diag(c((0.4+delta)**2,0.4**2)))


df<-expand.grid(seq(-3,3,len=100),seq(-3,3,len=100))
df_matrix<-data.matrix(df)
l<-lapply(seq_len(nrow(df_matrix)), function(i) df_matrix[i,])

result_mat=list()
Data<-c()
for (i in seq_len(nrow(rbind(z1,z2)))){
  mu<-rbind(z1,z2)[i,]
  result<-sapply(l,pmvnorm,lower=-Inf,mean=mu,sigma=I)
  Data<-rbind(Data,result)
}
inc<-seq(-3,3,len=100)[2]-seq(-3,3,len=100)[1]


distmat<-as.matrix(dist(Data,method = "manhattan")*inc)

depth_result<-depth_CPD(distmat,num_permut=num_permut ,c=0.1)

#### Graph CPD
E1 = mstree(dist(Data,method = "manhattan")*inc,ngmax = 5)
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

