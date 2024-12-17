TVPRR_cosso <-function(X,y,type,lambdavec=exp(linspace(-6,20,n=50)),lambda2=.1,olsprior=0,
                 sweigths=1,oosX=1,kfold=5,silent=1,alpha=.5,nnlasso=0,starting_betas=c(),oldlambda1=c(),
                 adaptive=0,aparam=0,tol=10^(-2),maxit=15,lambdabooster=1,homo.param=.75,
                 sv.param=.75){
  TYPE=type
  #to do : 
  #-update both procedures
  #-get rid of d's     DONE
  #-check for other compatibility issues with new versions
  #-run un it on this PC
  
  ######### Data pre-processing ###################
  
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
  if(is.null(starting_betas)==TRUE){
  if(length(lambdavec)>1){
    cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,lambdavec = lambdavec,k=kfold,plot=abs(silent-1),sweigths=sweigths)
    lambda1 <- cvlist$minimizer
    mse <- cvlist$minima
    
    #Final estimation
    grr<-dualGRR(ZZ,yy,dimX,lambda1,olsprior=olsprior
                 ,lambda2=lambda2,sweigths=sweigths)
    betas_grr <- grr$betas_grr
  }
  
  #CV
  if(length(lambdavec)==1){
    #Final estimation
    lambda1=lambdavec
    grr<-dualGRR(ZZ,yy,dimX,lambda1,lambda2=lambda2,sweigths=sweigths,olsprior = olsprior)
    betas_grr <- grr$betas_grr
  }
  }
  else{ #if a starting beta vector is given
    #need to re-normalize betas:
    coef = starting_betas
    for(j in 1:dim(X)[2]){
      coef[1,j+1,]=coef[1,j+1,]/scalingfactor[j]
    }
    coef[1,1,]=coef[1,1,]/sdy
    betas_grr=coef
    yhatold=yy
    for(tt in 1:dim(betas_grr)[3]){yhatold[tt]=betas_grr[1,,tt]%*%append(1,X[tt,])}
    grr=list(yhat=yhatold)
    lambda1 = oldlambda1
    }
  ######### TYPE 2 (TWO STEPS ARR) ############
  if(TYPE==2){
    umat = betas_grr[1,,2:nrow(ZZ)]-betas_grr[1,,1:(nrow(ZZ)-1)]
    sigmasq= diag(umat%*%t(umat))^homo.param
    sigmasq=(sigmasq/mean(sigmasq))
    
    e = yy - grr$yhat
    EW = (hush(garchFit(data=e,~ garch(1,1)))@sigma.t)^sv.param
    EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
    
    if(length(lambdavec)>1){cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,lambda2 = lambda2,eweigths = EW,
                            lambdavec = lambdavec,k=kfold,plot=abs(silent-1),sweigths=sigmasq)
    usethislambda1=cvlist$minimizer
    }
    else{usethislambda1=lambdavec}
    
    grrats<-dualGRR(ZZ,yy,dimX,lambda1=usethislambda1,eweigths = EW,olsprior = olsprior,
                      lambda2=lambda2,sweigths = sigmasq)
    betas_grrats <- grrats$betas_grr

    grra <- list() #since we did not do it
  }
  
  ######### TYPE 3 (ITERATIVE COSSO) ############
  if(TYPE==3){
    grrats <- list()
    maxit=maxit
    tol=tol

    #SV
    e = yy - grr$yhat
    e[(e^2)>10*sd(e)]=0 #glicth detector
    garchstuff = hush(garchFit(data=e,~ garch(1,1)))
    EW = (garchstuff@sigma.t)^sv.param
    EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
    
    msevec <- c() #mse1
    sigmavec <-c()
    diffvec <- c()

    betas_grra = betas_grr
    lambda1pos = which(lambdavec==lambda1)
    lambda1vec_subset = max(1,lambda1pos-3):min(length(lambdavec),lambda1pos+3)
    
    step=1
    par(mfrow=c(1,1))
    gamma= lambda1
    repeat {
      
      #Update weigths
      gammareg= gamma_reg(yy,X,betas_grra[1,,],adaptive=adaptive,aparam=aparam,
                          silent=silent,lambdabooster=lambdabooster,eweigths=EW)
      gamma_new=gammareg$gamma
      gamma_new = gamma_new[(dimX+1):(dimX*2)] 
      #print(gammareg$mse)
      diff = mean(abs(gamma_new-gamma))
      diffvec <- append(diffvec,diff)
      gamma=gamma_new
      
      #LASSO DOES STANDIZATION???
      if(silent==0){barplot(t(gamma))}
      
      #Update paths
      cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,lambda2 = lambda2,eweigths = EW,
                             lambdavec = lambdavec[lambda1vec_subset],
                             k=kfold,plot=abs(silent-1),sweigths=gamma)
      
      grra<-dualGRR(ZZ,yy,dimX,lambda1=cvlist$minimizer,eweigths=EW,
                          lambda2=lambda2,sweigths = gamma)
      
      #SV update
      e = yy - grra$yhat
      e[(e^2)>10*sd(e)]=0 #glicth detector
      garchstuff = hush(garchFit(data=e,~ garch(1,1)))
      EW = (garchstuff@sigma.t)^sv.param
      EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
      
      betas_grra=grra$betas_grr
      
      step=step+1

      if(diff<tol){
        print(paste('converged in', step,'steps'))
        break
      }
      if(sum(gamma)==0){
        print(paste('converged in', step,'steps because all variances are 0'))
        break
      }
      if(step==maxit){
        print('no convergence')
        break
      }
    }
    

    ##DO A LAST GAMMA STEP TO DO LASSO ON B_0

}
  fcast <- c()
  if(TYPE==1){coef <- grr$betas_grr}
  if(TYPE==2){coef <- grrats$betas_grr}
  if(2<TYPE){coef <- grra$betas_grr}
  
  #re-scale 
  for(j in 1:dim(X)[2]){
    coef[1,j+1,]=coef[1,j+1,]*scalingfactor[j]
  }
  coef[1,1,]=coef[1,1,]*sdy
  
  if(TYPE==1){grr$betas_grr<-coef
  grr$yhat =  sdy*grr$yhat}
  if(TYPE==2){grrats$betas_grr<-coef
  grrats$yhat =  sdy*grrats$yhat}
  if(2<TYPE){grra$betas_grr<-coef
  grra$yhat =  sdy*grra$yhat}
  
  if(silent==0){
    if(TYPE==3){plot.ts(cbind(msevec,append(diffvec,c(0,0))))}
    
    #what does beta look like
    par(mfrow=c(3,3))
    for(j in 1:dimX){ts.plot(betas_grr[1,j,])}
    par(mfrow=c(1,1))
    par(mfrow=c(3,3))
    for(j in 1:dimX){ts.plot(coef[1,j,])}
    par(mfrow=c(1,1))
  }
  
  if(length(oosX)>1){
    fcast <- t(as.matrix(append(1,oosX))) %*% as.matrix(coef[1,,dim(coef)[3]]) #basically the last beta_t
  }
  
  return(list(grr=grr,grrats=grrats,grra=grra,fcast=fcast,mse=mse,lambda1=lambda1))
}


gamma_reg = function(y,X,betas,adaptive=0,aparam=0,silent=1,lambdabooster=1,eweigths=1){
  Xplus = cbind(1,X)
  Zmod= t(betas[,])*Xplus
  if(adaptive==1){Zmod=reweighter(Zmod,betas,aparam=aparam)}
  Zall = cbind(Xplus,Zmod)
  w=rep(1,dim(Zall)[2])
  w[1:dim(Xplus)[2]]=0 #1/dim(Zall)[1]
  CV_nnlasso = cv.glmnet(y=y,x=Zall,lower.limits=0,intercept=0,
                         penalty.factor=w,standardize=FALSE,
                         lambda=exp(seq(log(0.001), log(5), length.out=100)),weights=eweigths)
  if(silent==0){plot(CV_nnlasso)}
  mse = min(CV_nnlasso$cvm)
  ll = rep(0,dim(Zall)[2])
  ll[1:dim(Xplus)[2]]=-10000
  fit = glmnet(y=y,x=Zall,lower.limits=ll,lambda=CV_nnlasso$lambda.min*lambdabooster,standardize=FALSE,
               penalty.factor=w,intercept=0,weights=eweigths)
  gamma = as.vector(fit$beta)
return(list(mse=mse,gamma=gamma))
}

reweighter = function(X,betas,aparam=0){ 
  umat = betas[,2:dim(betas)[2]]-betas[,1:(dim(betas)[2]-1)]
  sigmasq= diag(umat%*%t(umat))/dim(betas)[2]
  w=(sigmasq)^aparam #/mean(sigmasq)
  w[w==Inf]=0
  newX=X

  for(m in 1:dim(X)[2]){
    newX[,m]= w[m]*X[,m]
  }
  return(newX)
}

#to silence the GARCH function
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

#101
#to get back 100