TVPRR <-function(X,y,type,lambdavec=exp(linspace(-6,20,n=50)),lambda2=.1,olsprior=0,opt.sige=1,
                 sweigths=1,oosX=1,kfold=5,silent=1,alpha=.5,nnlasso=0,tol=10^-5,maxit=20){
  TYPE=type
  #to do : 
  #-update both procedures
  #-get rid of d's     DONE
  #-check for other compatibility issues with new versions
  #-run un it on this PC
  
  mse<-c() #
  grrats<-list()
  grra<- list()
  fcast<- c()
  
  #Standardize variances
  scalingfactor = rep(1,dim(X)[2])
  
  for(j in 1:dim(X)[2]){
    scalingfactor[j]=sd(y)/sd(X[,j])
    X[,j]=(X[,j])/sd(X[,j])
  }
  sdy =sd(y)
  y = (y)/sd(y)
  
  ZZ <- Zfun(X) #Basis expansion
  yy <- y #univariate model
  dimX= dim(X)[2]+1 # + constant
  
  ######### TYPE 1 (FOR THE BASE MODEL OR ESTIMATION GIVEN HPS) ############
  #CV
  if(length(lambdavec)>1){
  cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,lambdavec = lambdavec,k=kfold,plot=abs(silent-1),sweigths=sweigths)
  lambda1 <- cvlist$minimizer
  mse <- cvlist$minima
  
  #Final estimation
  grr<-dualGRR(ZZ,yy,dimX,lambda1,lambda2=lambda2,sweigths=sweigths,olsprior = olsprior)
  betas_grr <- grr$betas_grr
  }
  
  #CV
  if(length(lambdavec)==1){
    #Final estimation
    lambda1=lambdavec
    grr<-dualGRR(ZZ,yy,dimX,lambda1,lambda2=lambda2,sweigths=sweigths,olsprior = olsprior,eweigths = opt.sige)
    betas_grr <- grr$betas_grr
    }

  oldbetas_grr = betas_grr
  ######### TYPE 2 (TWO STEPS ARR) ############
  if(TYPE==2){
    umat = betas_grr[1,,2:nrow(ZZ)]-betas_grr[1,,1:(nrow(ZZ)-1)]
    sigmasq= diag(umat%*%t(umat))
    sigmasq=sigmasq/mean(sigmasq)
    
    #cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,lambda2 = lambda2,
    #                        lambdavec = lambda1,k=kfold,plot=abs(silent-1),sweigths=sigmasq)
    #lambda1 <- cvlist$minimizer
    #mse <- cvlist$minima
    msev=c()
    newlambdavec = linspace(0.1*lambda1,lambda1,n=10)
    for(j in 1:length(newlambdavec)){
    grrats<-dualGRR(ZZ,yy,dimX,newlambdavec[j],lambda2=lambda2,sweigths = sigmasq)
    betas_grrats <- grrats$betas_grr
    
    if(nnlasso==1){
      #Generate regressors
      Xplus = cbind(1,X)
      Zmod= t(betas_grrats[1,,])*Xplus
      PX=Xplus%*%solve(crossprod(Xplus)+lambda2*diag(dimX))%*%t(Xplus)
      yminX = (diag(length(y))-PX)%*%y
      MXZmod = (diag(length(y))-PX)%*%Zmod
      Zall = cbind(Xplus,Zmod)
      CV_nnlasso = cv.glmnet(y=y,x=Zall,lower.limits=0,intercept=0)
      mse = min(CV_nnlasso$cvm)
      ll = rep(0,dim(Zall)[2])
      ll[1:dimX]=-10000
      gamma = as.matrix(glmnet(y=y,x=Zall,lower.limits=ll,lambda=CV_nnlasso$lambda.min,
                               intercept=0)$beta)
      print(mse)
      msev = cbind(msev,mse)
      if(silent==0){
        plot(CV_nnlasso)
      }
    }    
    }
    #mse = min(msev)


    grra <- list()
  }
  
  ######### TYPE 3 (ITERATIVE ARR) ############
  if(2<TYPE & TYPE<6){
    grrats <- list()
    maxit=maxit
    tol=tol
    gamma=2
    q=1.45 #is actually doing l_0 in this paper's case
    eweigths=1 #for now
    
    lbeta=dim(betas_grr)[3]
    umat = betas_grr[1,,2:lbeta]-betas_grr[1,,1:(lbeta-1)]
    sigmasq= diag(umat%*%t(umat)) #/nrow(umat) #/lbeta
    sigmasq_ts = sigmasq
    if(silent==0){barplot(sigmasq)}
    sigmasqor <- sigmasq
    alpha=alpha
    
    msevec <- c() #mse1
    sigmavec <-c()
    diffvec <- c()
    gvcond <-mean(sigmasq)
    invnewlambda=sigmasq
    
    step=1
    par(mfrow=c(1,1))
    repeat {
      sigmavec <- rbind(sigmavec,sigmasq)
      #stack <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,lambda2=lambda2,
      #                      lambdavec = lambda1,k=3,plot=1,sweigths=invnewlambda/gvcond,eweigths=eweigths)
      #mse <-stack$minima
      #msevec <- append(msevec,mse)
      
      betas_grra<-dualGRR(ZZ,yy,dimX,lambda1=lambda1,olsprior = olsprior,
                          lambda2=lambda2,sweigths = invnewlambda/gvcond,eweigths = eweigths)$betas_grr
      step=step+1
      
      umat = betas_grra[1,,2:lbeta]-betas_grra[1,,1:(lbeta-1)]
      delta = 10^(-5*sqrt(step)) #(-2*sqrt(step))
      sigor <- diag(umat%*%t(umat))
      sigmasqnew = ((diag(umat%*%t(umat)))^gamma+delta^gamma)^((2-q)/gamma) #num stab /lbeta
      sigmasqnew[(sigor/max(sigor))<0.05]=0
      if(silent==0){barplot(sigor)
      print(mean(sigor))}
      diff = mean(abs(sigmasqnew-sigmasq))
      diffvec <- append(diffvec,diff)
      sigmasq=sigmasqnew

      newlambda = alpha*(1/sigmasq_ts)+(1-alpha)*(1/sigmasqnew) #sigmasq_ts
      invnewlambda = 1/newlambda
      gvcond <- mean(sigmasq) #(mean(sigmasq)/mean(sigmasqnew))^-1 #mean(sigmasq)
      
      if(diff<tol){
        print(paste('converged in', step,'steps'))
        break
      }
      if(step==maxit){
        print('no convergence')
        break
      }
    }
    
    sigmasqar = sigmasq
    
    #sigmasqstar <-sigmavec[which.min(msevec),]
    #print(paste('best step is',which.min(msevec), 'at mse',min(msevec)))
    
    if(silent==0){
    plot.ts(cbind(msevec,append(diffvec,c(0,0))))
    
    #what does beta look like
    par(mfrow=c(3,3))
    for(j in 1:dimX){ts.plot(betas_grr[1,j,])}
    par(mfrow=c(1,1))
    par(mfrow=c(3,3))
    for(j in 1:dimX){ts.plot(betas_grra[1,j,])}
    par(mfrow=c(1,1))
    }
    
    #trimmer of the best sigmasq
    if(TYPE==4){
    msevec <- c()
    sigmavec <-c()
    
    for(i in 1:dimX){
      sigmasqstari <-sigmasqstar
      sigmasqstari[order(sigmasq)<i]=0
      print(sigmasqstari)
      sigmavec <- rbind(sigmavec,sigmasqstari)
      
      stack <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,lambda2=lambda2,
                            lambdavec = lambda1,k=kfold,plot=1,sweigths=sigmasqstari,eweigths=eweigths)
      mse <-stack$minima
      msevec <- append(msevec,mse)
    }
    
    sigmasqstar <-sigmavec[which.min(msevec),]
    if(silent==0){
    print(paste('best step is',which.min(msevec), 'at mse',min(msevec)))
    print(sigmasqstar)
    ts.plot(msevec)
    }
    
    #final
    if(sum(sigmasqstar)!=0){SW=sigmasqstar/mean(sigmasqstar)}
    else{SW=0}
    
    cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,
                           lambdavec = lambda1,lambda2=lambda2,
                           k=kfold,plot=1,sweigths=SW)
    
    lambda1 <- cvlist$minimizer
    mse <- cvlist$minima
    
   
    grra<-dualGRR(ZZ,yy,dimX,lambda1,olsprior = olsprior,
                  lambda2=lambda2,sweigths = SW)
    
    }
    else{
      if(sum(sigmasqar)!=0){SW=sigmasqar/mean(sigmasqar)}
      else{SW=0}
      
      cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,
                             lambdavec = lambda1,lambda2=lambda2,
                             k=kfold,plot=1,sweigths=SW)
      
      lambda1 <- cvlist$minimizer
      mse <- cvlist$minima
      
      grra<-dualGRR(ZZ,yy,dimX,lambda1,olsprior=olsprior,
                    lambda2=lambda2,sweigths = SW)
    
      }
  }
  
  ######### TYPE 6 FACTOR TVP ############
  if(TYPE==6){
    print('spinning to factor space')
    grrats <- list()
    maxit=20
    tol=10^(-4)
    gamma=2
    q=1.45 #is actually doing l_0 in this paper's case
    eweigths=1 #for now
    
    lbeta=dim(betas_grr)[3]
    umat = betas_grr[1,,2:lbeta]-betas_grr[1,,1:(lbeta-1)]
    sigmasq= diag(umat%*%t(umat)) #/nrow(umat) #/lbeta
    sigmasq_ts = sigmasq
    if(silent==0){barplot(sigmasq)}
    sigmasqor <- sigmasq
    alpha=alpha
    
    msevec <- c() #mse1
    sigmavec <-c()
    diffvec <- c()
    gvcond <-mean(sigmasq)
    invnewlambda=sigmasq
    
    aa = fastZrot(oldbetas_grr,X)

    step=1
    par(mfrow=c(1,1))
    repeat {
      sigmavec <- rbind(sigmavec,sigmasq)
      
      betas_grra<-dualGRR(aa$Zrot,yy,dimX,lambda1=lambda1,
                          lambda2=lambda2,sweigths = invnewlambda/gvcond,eweigths = eweigths)$betas_grr
      step=step+1
      
      umat = betas_grra[1,,2:lbeta]-betas_grra[1,,1:(lbeta-1)]
      delta = 10^(-5*sqrt(step)) #(-2*sqrt(step))
      sigor <- diag(umat%*%t(umat))
      sigmasqnew = ((diag(umat%*%t(umat)))^gamma+delta^gamma)^((2-q)/gamma) #num stab /lbeta
      sigmasqnew[(sigor/max(sigor))<0.05]=0
      if(silent==0){barplot(sigor)
        print(mean(sigor))}
      diff = mean(abs(sigmasqnew-sigmasq))
      diffvec <- append(diffvec,diff)
      sigmasq=sigmasqnew
      
      newlambda = alpha*(1/sigmasq_ts)+(1-alpha)*(1/sigmasqnew) #sigmasq_ts
      invnewlambda = 1/newlambda
      gvcond <- mean(sigmasq) #(mean(sigmasq)/mean(sigmasqnew))^-1 #mean(sigmasq)
      
      if(diff<tol){
        print(paste('converged in', step,'steps'))
        break
      }
      if(step==maxit){
        print('no convergence')
        break
      }
    }
    
    sigmasqar = sigmasq
    
    #NEED A RECONSTRUCTION STEP OR FUNCTION
    
    
    #sigmasqstar <-sigmavec[which.min(msevec),]
    #print(paste('best step is',which.min(msevec), 'at mse',min(msevec)))
    
    if(silent==0){
      plot.ts(cbind(msevec,append(diffvec,c(0,0))))
      
      #what does beta look like
      par(mfrow=c(3,3))
      for(j in 1:dimX){ts.plot(betas_grr[1,j,])}
      par(mfrow=c(1,1))
      par(mfrow=c(3,3))
      for(j in 1:dimX){ts.plot(betas_grra[1,j,])}
      par(mfrow=c(1,1))
    }
    
    if(sum(sigmasqar)!=0){SW=sigmasqar/mean(sigmasqar)}
    else{SW=0}
    
    cvlist <- cvgs.bhk2015(ZZ=aa$Zrot,Y=yy,dimX = dimX,
                           lambdavec = lambda1,lambda2=lambda2,
                           k=kfold,plot=1,sweigths=SW)
    
    lambda1 <- cvlist$minimizer
    mse <- cvlist$minima
    
    grra<-dualGRR(aa$Zrot,yy,dimX,lambda1,
                  lambda2=lambda2,sweigths = SW)
    
    grra$betas_grr = FtoB(betas=grra$betas_grr,Lambda=aa$Lambda)$betas
  
    }

  
  
  fcast <- c()
  
  if(TYPE==1){coef <- grr$betas_grr}
  if(TYPE==2){coef <- grrats$betas_grr}
  if(2<TYPE){coef <- grra$betas_grr}
  
  #poor man's factor prior (18/10/09)
  #if(TYPE==1){coef <-RRankU(grr$betas_grr,X,y)$newbetas}
  #if(TYPE==2){coef <- RRankU(grrats$betas_grr,X,y)$newbetas}
  #if(2<TYPE){coef <- RRankU(grra$betas_grr,X,y)$newbetas}
  
  #print(dim(grra$betas_grr))
  #re-scale 
  for(j in 1:dim(X)[2]){
    coef[1,j+1,]=coef[1,j+1,]*scalingfactor[j]
  }
  coef[1,1,]=coef[1,1,]*sdy
  
  if(TYPE==1){grr$betas_grr<-coef}
  if(TYPE==2){grrats$betas_grr<-coef}
  if(2<TYPE){grra$betas_grr<-coef}
  
  if(length(oosX)>1){
    fcast <- t(as.matrix(append(1,oosX))) %*% as.matrix(coef[1,,dim(coef)[3]]) #basically the last beta_t
  }
  
  return(list(grr=grr,grrats=grrats,grra=grra,fcast=fcast,mse=mse,lambda1=lambda1))
}
    
RRankU <- function(betas,XX,yy){
  lbeta=dim(betas)[3]
  t=lbeta
  
  umat = betas[1,,2:lbeta]-betas[1,,1:(lbeta-1)]
  omega = umat%*%t(umat)
  ev<- eigen(omega)$values/sum(eigen(omega)$values)
  cev = cumsum(ev)
  
  p= dim(umat)[1]
  fanal = factor(t(umat),dim(umat)[1]) #factor analysis 
  newev = array(0,dim=c(1,p)) 
  evnotzero = eigen(omega)$values[(cev)<0.85] #selection rule for now
  if(length(evnotzero)>0){newev[1,1:length(evnotzero)] =  evnotzero} #put zeros on the other eigenvalues
  else{newev[1,1] =  eigen(omega)$values[1]}
  
  newLambda = array(0,dim=c(p,p))
  newLambda[,1:length(evnotzero)] = fanal$lambda[,1:length(evnotzero)] #lambda with imposed restrictions
  #print(length(evnotzero))
  #print(cev)
  
  F = t(fanal$factors)
  newU = newLambda%*%F
  
  #RECREATE BETA
  #GENERATE NEW Y WITH A LOOP
  # RUN LASSO VERSION 
  
  dimX = dim(betas)[2]
  newbetas = betas #array(NA,dim=c(1,dimX,dim(newU)[2]))
  nsF = newbetas
  for(k in 1:dimX){ #variables
    for(tt in 1:(dim(newU)[2])){
      nsF[1,k,tt] = sum(F[k,1:tt])
    }
  }
  
  eq=1
  #Recover betas from u's
  for(k in 1:dimX){ #variables
    for(tt in 2:(dim(newU)[2])){
      newbetas[eq,k,tt] = betas[1,k,1]+sum(newU[k,1:(tt-1)])
    }
  }
  
  newyy=yy
  partialyhat = yy
  for(tho in 1:t){
    partialyhat[tho]=(append(1,XX[tho,]))%*%(as.vector(newbetas[1,,tho]))
    newyy[tho]=yy[tho]-partialyhat[tho]
  }
  
  return(list(newbetas=newbetas,newyy =newyy,nsF=nsF))
}

Zrot<- function(betas,Z,dimX){
  
  lbeta=dim(betas)[3]
  t=lbeta
  
  umat = betas[1,,2:lbeta]-betas[1,,1:(lbeta-1)]
  omega = umat%*%t(umat)
  ev<- eigen(omega)$values/sum(eigen(omega)$values)
  cev = cumsum(ev)

  p= dim(umat)[1]
  fanal = factor(t(umat),dim(umat)[1]) #factor analysis 
  newev = array(0,dim=c(1,p)) 
  evnotzero = eigen(omega)$values[(cev)<1.2] #(selection rule for now) no seledtion here
  if(length(evnotzero)>0){newev[1,1:length(evnotzero)] =  evnotzero} #put zeros on the other eigenvalues
  else{newev[1,1] =  eigen(omega)$values[1]}
  newLambda = array(0,dim=c(p,p))
  newLambda[,1:length(evnotzero)] = fanal$lambda[,1:length(evnotzero)] #lambda with imposed restrictions
  
  I_Tmin = diag(t-1)
  Zrot <- Z[,1:(dim(Z)[2]-dimX)]%*%(kron(newLambda,I_Tmin))
  return(list(Zrot=Zrot,Lambda=newLambda))
}

recoverer = function(F,Lambda,betas){

newU = Lambda%*%F[1,,]
dimX = dim(F)[2]
newbetas = betas

nsF = array(NA,dim=c(1,dimX,dim(newU)[2]))
for(k in 1:dimX){ #variables
  for(tt in 1:(dim(newU)[2])){
    nsF[1,k,tt] = sum(F[1,k,1:tt])
  }
}

#Recover betas from u's
for(k in 1:dimX){ #variables
  for(tt in 2:(dim(newU)[2])){
    newbetas[1,k,tt] = betas[1,k,1]+sum(newU[k,1:(tt)])
  }
}

return(list(nsF=nsF,betas=newbetas))
}

fastZrot<- function(betas,X){
  
  lbeta=dim(betas)[3]
  t=lbeta
  
  umat = betas[1,,2:lbeta]-betas[1,,1:(lbeta-1)]
  omega = umat%*%t(umat)
  ev<- eigen(omega)$values/sum(eigen(omega)$values)
  cev = cumsum(ev)
  
  p= dim(umat)[1]
  fanal = factor(t(umat),dim(umat)[1]) #factor analysis 
  print(mean(diag(omega)))
  
  print(mean(diag(t(fanal$factors)%*%fanal$factors)))
  newev = array(0,dim=c(1,p)) 
  evnotzero = eigen(omega)$values[(cev)<.98] #(selection rule for now) no seledtion here
  if(length(evnotzero)>0){newev[1,1:length(evnotzero)] =  evnotzero} #put zeros on the other eigenvalues
  else{newev[1,1] =  eigen(omega)$values[1]}
  
  newLambda = array(0,dim=c(p,p))
  newLambda[,1:length(evnotzero)] = fanal$lambda[,1:length(evnotzero)]/dim(betas)[2] #lambda with imposed restrictions
  #print(newLambda)
  #I_Tmin = diag(t-1)
  #Zrot <- cbind(Z[,1:(dim(Z)[2]-dimX)]%*%(kron(newLambda,I_Tmin)),Z[,(dim(Z)[2]-dimX+1):dim(Z)[2]])
  
  XX= cbind(matrix(1,(nrow(X)),1),X)
  Z=array(66,dim=c((nrow(X)),(nrow(X)-1),(ncol(X)+1)))
  
  #The sequence of tau's
  param_seq  <- seq(2,nrow(X),by=1)
  
  #New regressors
  for(tt in 1:(dim(X)[1])){
    for(tho in param_seq){
      if(tho<=tt){
        for(jj in 1:(ncol(X)+1)){
          Z[tt,tho-1,jj]= newLambda[,jj]%*%XX[tt,]
        }
      }
      else{
        Z[tt,tho-1,]=0
      }
    }
  }
  
  Zprime <- c()
  for(kk in 1:(ncol(X)+1)){
    Zprime = cbind(Zprime,Z[,,kk])
  }
  
  Zrot <- cbind(Zprime,XX)
  return(list(Zrot=Zrot,Lambda=newLambda))
}

FtoB <- function(Lambda,betas){
  
  F = betas[1,,2:dim(betas)[3]]-betas[1,,1:(dim(betas)[3]-1)]
  newU = Lambda%*%F
  
  #RECREATE BETA
  #GENERATE NEW Y WITH A LOOP
  # RUN LASSO VERSION 
  
  dimX = dim(betas)[2]
  newbetas = betas #array(NA,dim=c(1,dimX,dim(newU)[2]))
  nsF = newbetas
  nsF[1,,1]=F[,1]
  for(k in 1:dimX){ #variables
    for(tt in 2:(dim(betas)[3])){
      nsF[1,k,tt] = sum(F[k,1:(tt-1)])
    }
  }
  
  eq=1
  #Recover betas from u's
  for(k in 1:dimX){ #variables
    for(tt in 2:(dim(betas)[3])){
      newbetas[eq,k,tt] = betas[1,k,1]+sum(newU[k,1:(tt-1)])
    }
  }
  return(list(betas=newbetas,nsF=nsF))
}