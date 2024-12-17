#simul 
{    
  heta <- array(NA, dim=c(obs,M+1))
  X <- array(NA, dim=c(obs,M-1))
  
  for(j in 1:(M+1)){
    set.seed(1234+s*1000*setup+j)
    heta[,j]<- as.matrix(rnorm(obs))
  }
  for(j in 1:(M-1)){
    set.seed(1234+s*1000*setup+j*3)
    X[,j]<- as.matrix(rnorm(obs))/sqrt(M/20)
  }
  eps       <- as.matrix(rnorm(obs))
  y = matrix(0,obs,1)
  yhat=y
  beta = matrix(0,obs,M+1)
  
  for(j in 1:(M+1)){
    beta[1,j]=0 #don't change that for the sin()
  }
  
  Z=array(NA,dim=c(obs,M+1))
  Z[,1]=1
  howmanytvp=M+1
  truesig = array(0,dim=c(M+1,1))
  
  for(t in 2:obs) {
    for(j in 1:howmanytvp){
      if(j>0.5*fracvec[frac]*howmanytvp){beta[t,j]=compressor*(sin(t/20) + sig_heta*heta[t,j])
      if(j>fracvec[frac]*howmanytvp){beta[t,j]=compressor*(((-1)^j)*sqrt(t/(300)) + sig_heta*heta[t,j])
      truesig[j] = 1
      
      }} #(1/j)*
      #if((-1)^j>0){beta[t,j]=compressor*(sin(t/20) + sig_heta*heta[t,j])} #(1/j)*
    }
  }
  if(var(beta[,2])!=0){
    beta[,2]= (beta[,2]-min(beta[,2]))/(max((beta[,2]))-min(beta[,2]))
    beta[,2]= beta[,2]*0.4+0.2
  }
  for(t in 2:obs) {
    Z[t,2]= y[t-1]
    Z[t,3:(M+1)]=X[t-1,]
    yhat[t]=Z[t,]%*%(beta[t,])
    y[t] <-Z[t,]%*%(beta[t,]) +compressor*sig_e[t]*eps[t]
  }
  
  # Rename data
  usmacro <- cbind(y,X)
  #usmacro=usmacro[50:nrow(usmacro),]
  T=nrow(usmacro)
  beta=beta[2:nrow(beta),]
  beta = t(beta)
  
}