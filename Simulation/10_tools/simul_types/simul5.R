#simul 
{    
  heta <- array(NA, dim=c(obs,M+1))
  X <- array(NA, dim=c(obs,M-1))
  
  for(j in 1:(M+1)){
    set.seed(1234+56*s*1000*setup+j*s)
    heta[,j]<- as.matrix(rnorm(obs))
  }
  for(j in 1:(M-1)){
    set.seed(1234+s*1000*setup+j*3)
    X[,j]<- as.matrix(rnorm(obs))/sqrt(M/20)
  }
  
  epsall       <- as.matrix(rnorm(obs*2))
  eps =epsall[1:obs]
  eps2 =epsall[(obs+1):(2*obs)] #/2
  
  y = matrix(0,obs,1)
  beta = matrix(0,obs,M+1)
  
  for(j in 1:(M+1)){
    beta[1,j]=rnorm(1)/2 #don't change that for the sin()
  }
  
  Z=array(NA,dim=c(obs,M+1))
  Z[,1]=1
  howmanytvp=M+1
  truesig = array(0,dim=c(M+1,1))
  yhat=y 
  loadings=cbind(rnorm(M+1),rnorm(M+1),rnorm(M+1))
  jump =rep(0,obs)
  jump[200]=.4
  facs = t(cbind(0*jump+1*rnorm(obs,0,0.011),0.03*append(rep(-.12,round(obs/2)+1),rep(.12,obs/2)),0.015*sin(0.05*1:T)))
  yhat2=yhat
  y2=y
  
  for(t in 2:obs) {
    for(j in 1:howmanytvp){
      if(j>fracvec[frac]*howmanytvp){beta[t,j]=beta[t-1,j] + compressor*(loadings[j,]%*%facs[,t]) #+ sig_heta*heta[t,j]
      #if(j>fracvec[frac]*howmanytvp){beta[t,j]=compressor*(((-1)^j)*sqrt(t/(100)) + sig_heta*heta[t,j])
      truesig[j] = 1
      
      }
      else{beta[t,j]=beta[t-1,j]}
      } #(1/j)*
      #if((-1)^j>0){beta[t,j]=compressor*(sin(t/20) + sig_heta*heta[t,j])} #(1/j)*
    
    Z[t,2]= y[t-1]
    Z[t,3:(M+1)]=X[t-1,]
    yhat[t]=Z[t,]%*%(beta[t,])
    y[t] <-Z[t,]%*%(beta[t,]) +compressor*sig_e[t]*eps[t]
    
    Z[t,2]= y2[t-1]
    Z[t,3:(M+1)]=X[t-1,]
    yhat2[t]=-0.8*Z[t,]%*%(beta[t,])
    y2[t] <-yhat2[t] +compressor*sig_e[t]*eps2[t]
  }
  
  # Rename data
  usmacro <- cbind(y,X)
  #usmacro=usmacro[50:nrow(usmacro),]
  T=nrow(usmacro)
  beta=beta[2:nrow(beta),]
  beta = t(beta)
  
  #print(summary(lm(y~yhat)))
  
}