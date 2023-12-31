require(pracma)
for (p in c(30,90,180)){
  M<-(matrix(1,1,p))/sqrt(p)
  N<-nullspace(M)
  U<-cbind(t(M),N)
  
  v<-c()
  for (k in (1:p)){
    v<-c(v,cos(k*pi/p)+1.5)
  }
  Delta<-diag(v)
  Sigma<-(U) %*% Delta %*% t(U)
  
  #save(Sigma, file="Sigma_90.RData")
  path<-paste("Sigma_",p,'.Rdata',sep="")
  save(Sigma, file=path)
}
