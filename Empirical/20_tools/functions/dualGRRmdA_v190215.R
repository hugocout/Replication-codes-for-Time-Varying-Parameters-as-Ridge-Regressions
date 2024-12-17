dualGRR <-function(Zprime,y,dimX,lambda1,lambda2=0.001,CI=0,sweigths=1,olsprior=1,eweigths=1,GCV=0,calcul_beta=1,nf=0){
  ### Inputs are: ####
  #Zprime is the expanded X's
  #y is self-explanatory 
  #lambda[1] is the the ratio of sigma_u and sigma_epsilon
  #lambda[2] is the penalty for non-u's parameters
  # d is the derivative being penalized (1 or 2)
  #####################
  
  
  ############ PREPARE PRIOR MATRIX ############
  lambdas = c(lambda1,lambda2)
  if(nf==0){nf=dimX}
  
  #lambda[2]  cannot be too small. let's impose a safety guard.
  if(lambdas[2]<0.0000001){
    print('Lambda_2 imposed to be at least 0.0001')
    lambdas[2]<-0.0000001
  }
  
  #create inv(Lambda_p) implicitely 
  ncolZ = ncol(Zprime)
  T=(dim(Zprime)[2]-dimX)/nf +1

  if(length(sweigths)==1){sweigths=rep(1,1,nf)}
    Kmat_half<-t(Zprime)
    for(m in 1:(nf)){
      begin=(m-1)*(T-1)+1
      end=(m)*(T-1)
      Kmat_half[begin:end,]= (1/lambdas[1])*sweigths[m]*Kmat_half[begin:end,]
      if(-nf>-dimX){Kmat_half[begin,]=1000*Kmat_half[begin,]}
    }
  Kmat_half[(ncolZ-dimX+1):ncolZ,]=(1/lambdas[2])*Kmat_half[(ncolZ-dimX+1):ncolZ,]
  
  ############ DUAL GRR ############
  
  Lambda_T<-diag(nrow(Zprime))
  if(length(eweigths)>1){Lambda_T<- diag(as.numeric(eweigths))}
  
  param=nf*(T-1) + dimX
  obs = nrow(Zprime)
  
  if(param>obs){ #Go Dual
  Kmat_du <- Zprime%*%Kmat_half
  
  if(olsprior==0){
    alpha <- solve(Kmat_du + Lambda_T,y)
    uhat <- Kmat_half%*%alpha
  }
  if(olsprior==1){
    X= Zprime[,(ncolZ-dimX+1):ncolZ]
    beta_ols<-solve(crossprod(X),crossprod(X,y))
    alpha <- solve(Kmat_du + Lambda_T,(y-X%*%beta_ols))
    uhat <- Kmat_half%*%alpha
    uhat[(ncolZ-dimX+1):ncolZ]= uhat[(ncolZ-dimX+1):ncolZ]+beta_ols
  }
  yhat <- Kmat_du%*%alpha
  }
  
  else{ #Go Primal
    for(tt in 1:obs){Kmat_half[,tt]=Kmat_half[,tt]*eweigths[tt]^-1}
    Kmat_pri <- Kmat_half%*%Zprime
    
    if(olsprior==0){
      uhat <- solve(Kmat_pri + diag(param),Kmat_half%*%y)
    }
    if(olsprior==1){
      X= Zprime[,(ncolZ-dimX+1):ncolZ]
      beta_ols<-solve(crossprod(X),crossprod(X,y))
      uhat <- solve(Kmat_pri + diag(param),Kmat_pri%*%(y-X%*%beta_ols))
      uhat[(ncolZ-dimX+1):ncolZ]= uhat[(ncolZ-dimX+1):ncolZ]+beta_ols
    }
    yhat <- Zprime%*%uhat
  }
  
  ############ RECOVER BETA'S ############
  betas_grr <- array(0,dim=c(dim(as.matrix(y))[2],nf,T))
  
  if(calcul_beta==1){

    for(eq in 1:dim(as.matrix(y))[2]){
      for(k in 1:nf){betas_grr[eq,k,1] = uhat[(dim(uhat)[1]-dimX+k),eq]}
        for(t in 2:(T)){
          begin=(k-1)*(T-1)+1
          positions=c()
          for(k in 1:nf){positions = append(positions,(k-1)*(T-1)+(t-1))}
          #print(end)
          #if(d==1){begin=end}
          betas_grr[eq,,t] = betas_grr[eq,,t-1]+uhat[positions,eq]
        }
      }
  }
  ############ CI'S ############
if(CI>0){
  
  Kmh = Kmat_half
  if(param>obs){ #Go Dual
    for(tt in 1:obs){Kmh[,tt]=Kmat_half[,tt]*eweigths[tt]^-1}
    #Kmat_pri <- Kmh%*%Zprime
  }
  
  if(dimX>nf){
    Kmh=Kmh[1:(dim(Kmat_half)[1]-dimX),]
    #Kmat_pri = Kmat_pri[1:(dim(Kmat_half)[1]-dimX),1:(dim(Kmat_half)[1]-dimX)]
  }
  else{
  
  #Putting back Kmat_half into order
  #Put back the first col of each
  for(k in 1:dimX){
  Kmh[1+T*(k-1),] = Kmh[(dim(Kmat_half)[1]-dimX+k),]
  begin = 2+T*(k-1)
  end = T+T*(k-1)
  Kmh[begin:end,]= Kmh[(begin-k):(end-k),]
  }}
  #Getting #Vb
  
  invCov = chol2inv(chol(tcrossprod(Kmh)+diag(dim(Kmh)[1])))

  CT= array(1,dim=c(T-1,T-1))
  CT= lower.triangle(CT)
  C=kronecker(diag(nf),CT)
  #print(dim(invCov))
  #print(dim(C))
  


  #Vb=diag(C%*%tcrossprod(invCov,C))*sd(y-yhat)/sqrt(dim(Kmh)[2])
  #sandwich = invCov%*%Kmat_pri%*%invCov
  sandwich=invCov*(sd(y-yhat)^2)*(lambdas[1])^(-1)
  Vb=sqrt(diag(C%*%tcrossprod(sandwich,C))) #/(dim(Kmh)[2])
  
  betas_grr_ci <- array(0,dim=c(dim(as.matrix(y))[2],nf,T,2))
  for(c in 1:2){
    #temp<- uhat + ((-1)^c)*1.96*Vu
      for(eq in 1:dim(as.matrix(y))[2]){
        for(k in 1:nf){ #variables
          #betas_grr_ci[eq,k,1,c] = temp[(dim(uhat)[1]-dimX+k),eq]
          for(t in 1:(T)){
            #begin=(k-1)*(T-1)+1
            #end =(k-1)*(T-1)+(t-1)
            betas_grr_ci[eq,k,t,c] = betas_grr[eq,k,t]+((-1)^c)*CI*Vb[(k-1)*(T-1)+1]
          }
        }
    }
  }

  return(list(uhat=uhat,betas_grr=betas_grr,yhat=yhat,betas_grr_ci=betas_grr_ci,lambdas=lambdas,
                sigmasq=sweigths))
}
  else{return(list(uhat=uhat,betas_grr=betas_grr,yhat=yhat,sigmasq=sweigths,lambdas=lambdas))
}
}
