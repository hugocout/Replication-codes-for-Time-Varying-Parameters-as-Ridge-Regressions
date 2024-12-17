TVPRR_cossoF <-function(X,y,type,lambdavec=exp(linspace(-6,20,n=50)),lambda2=.1,
                 sweigths=1,oosX=1,kfold=5,silent=1,alpha=.5,nnlasso=0,olsprior=0,
                 adaptive=0,aparam=0,tol=10^(-2),maxit=15,starting_betas=c(),
                 lambdabooster=1,var.share=.8,override=0,id=3,homo.param=.75,
                 sv.param=.75,fp.model=1,max.step.cv=1){
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
  
  if(is.null(starting_betas)==TRUE){
  
  ######### TYPE 1 (FOR THE BASE MODEL OR ESTIMATION GIVEN HPS) ############
  #CV
  if(length(lambdavec)>1){
    cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=yy,dimX = dimX,lambdavec = lambdavec,k=kfold,plot=abs(silent-1),sweigths=sweigths)
    lambda1 <- cvlist$minimizer
    mse <- cvlist$minima
    
    #Final estimation
    grr<-dualGRR(ZZ,yy,dimX,lambda1,olsprior = olsprior
                 ,lambda2=lambda2,sweigths=sweigths)
    betas_grr <- grr$betas_grr
  }
  
  #CV
  if(length(lambdavec)==1){
    #Final estimation
    lambda1=lambdavec
    grr<-dualGRR(ZZ,yy,dimX,lambda1,lambda2=lambda2,sweigths=sweigths,olsprior=olsprior)
    betas_grr <- grr$betas_grr
  }
  lambda1_step1=lambda1
  }
  
  else{
    lambda1_step1=NA
    coef = starting_betas
    for(j in 1:dim(X)[2]){coef[1,j+1,]=coef[1,j+1,]/scalingfactor[j]}
    coef[1,1,]=coef[1,1,]/sdy
    betas_grr=coef
    yhatold=yy
    for(tt in 1:dim(betas_grr)[3]){yhatold[tt]=betas_grr[1,,tt]%*%append(1,(X[tt,]))}
    #for(tt in 1:length(yy)){yhatold[tt]=betas_grr[1,,tt]%*%append(1,(X[tt,]))}
    grr=list(yhat=yhatold)
  }
###### FACTOR COSSO #####
  if(TYPE==2){
    grrats <- list()
    maxit=maxit
    tol=tol

    #SV
    e = yy - grr$yhat
    e[(e^2)>10*mad(e)]=0 #glicth detector
    #print(e)
    garchstuff = hush(garchFit(data=e,~ garch(1,1)))
    EW = (garchstuff@sigma.t)^sv.param
    EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
    
    msevec <- c() #mse1
    sigmavec <-c()
    diffvec <- c()

    betas_grra = betas_grr
    
    step=1
    par(mfrow=c(1,1))

    bestmse=Inf
    for(nf in 1:length(override)){
    innerbestmse=Inf
    
    aa = BtoF(betas_grr[1,,],var.share=var.share,override = override[nf],id=id,vizualizar = abs(silent-1))
    newF =aa$f
    #print(aa$Lambda)
    newLambda = aa$Lambda
    newU=newLambda%*%newF
    gamma = 1
    SW=1
    
    repeat {
      
      r=dim(newLambda)[2]
      if(is.null(r)==TRUE){r=1}
      
      if(silent==0){
        #what does factor look like
        
        if(r==1){ts.plot(cumsum(newF))}
        else{
          par(mfrow=c(3,3))
          for(j in 1:dim(as.matrix(newLambda))[2]){ts.plot(cumsum(as.matrix(newF)[j,]))}
          par(mfrow=c(1,1))   
        }
      }

      #Update paths
      newZ=R_Lambda(X=X,newLambda = newLambda) #,Z=ZZ)
      
      if(-step>=-max.step.cv){
      stuff <- CV.KF.MV(ZZ=newZ,Y=yy,dimX = dimX,lambdavec = 
                                lambdavec,k=kfold,plot=abs(silent-1),eweigths=EW,
                              sweigths=SW,nf=dim(newLambda)[2],lambda2=lambda2)
      #if(step>1){lambda1=mean(append(stuff$minimizer,lambda1))}
      #if(step>1){lambda1=lambda1}
      #else{lambda1=stuff$minimizer}
      lambda1=stuff$minimizer
      msevec = append(msevec,stuff$minima)
      }
      
      grra<-dualGRR(Zprime=newZ,y=yy,dimX=dimX,lambda1=lambda1,eweigths=EW,
                    lambda2=lambda2,sweigths = SW,CI=0,olsprior=olsprior,
                    nf=dim(newLambda)[2])
      newF=grra$betas_grr[1,,2:dim(grra$betas_grr)[3]]-grra$betas_grr[1,,1:(dim(grra$betas_grr)[3]-1)]
      
      #SV update
      e = yy - grra$yhat
      #plot.ts(yy)
      #plot.ts(grra$yhat)
      e[(e^2)>10*mad(e)]=0 #glicth detector
     ##plot.ts(e)
      garchstuff = hush(garchFit(data=e,~ garch(1,1)))
      EW = (garchstuff@sigma.t)^sv.param
      EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
      
      #Normalization
      if(id==1){
      #Normalization of F and creation of newZ for loadings regression
      utilT=dim(grra$betas_grr)[3]-1
      if(length(newF)>utilT){
      orthomat = (tcrossprod(newF))
      SW = (diag(orthomat)/mean(diag(orthomat)))^(.5)
      orthoF= (pracma::sqrtm(orthomat)$Binv)%*%newF
      }
      else{orthoF=newF/sd(newF)}
      newF=orthoF
      }
      newZ = R_f2(f=newF,dimX=dimX,Z=ZZ)

      #Loadings regression
      gammareg= gamma_reg_f(y=yy,X=X,Z=newZ,adaptive=adaptive,aparam=aparam,eweigths=EW,alpha=alpha,
                          silent=silent,lambdabooster=lambdabooster,id=id)
      gamma_new=gammareg$gamma
      beta_0 = gamma_new[1:dimX]
      gamma_new = gamma_new[(dimX+1):(length(gamma_new))] 
      gamma=gamma_new
      #print(gammareg$Lambda)
      newLambda = gammareg$Lambda

      if(r!=1){
      dumpvec=c()
      for(rr in 1:r){if(-sum(abs(newLambda[,rr]))>-0.01*tol){dumpvec=append(dumpvec,rr)}}
      if(is.null(dumpvec)==FALSE){
        newLambda=as.matrix(newLambda[,-dumpvec])
        newF=as.matrix(newF[-dumpvec,])
        if(dim(newF)[1]>dim(newF)[2]){newF=t(newF)}
        r = r - length(dumpvec)
      }}
      else{if(-sum(abs(newLambda))>-0.01*tol){r=0}}

      if(silent==0){barplot(t(gamma))}

      #newLambda[newLambda==min(newLambda)]=0
      
      #Impose normalization of Lambda matrix
      oldU=newU
      newU=newLambda%*%newF
      
      if(r==0){
        print(paste('converged in', step,'steps because all variances are 0'))
        break
      }
      
      if(id==1){
      aa = UtoF((newU),var.share=var.share/2,override = override[nf],id=id,vizualizar = abs(silent-1))
      newF =aa$f
      newLambda = aa$Lambda
      if(r!=dim(newLambda)[2]){SW=1}
      }

      #Convergence crit check
      diff = mean(abs(newU-oldU))
      diffvec <- append(diffvec,diff)
      step=step+1
      
      #keeping the best mse model
      if(innerbestmse>tail(msevec,1)){
        innerbestmse=tail(msevec,1)
        inner.model.vault = list(newU=newU,beta_0=beta_0)
      }

      if(diff<tol){
        print(paste('converged in', step,'steps'))
        break
      }
      if(step==maxit){
        print('no convergence')
        break
      }
    }
    step=1

    #the FP or the best mse?
    if(fp.model==1){innerfinalmse=tail(msevec,1)}
    else{
      innerfinalmse=innerbestmse
      newU=inner.model.vault$newU
      beta_0 = inner.model.vault$beta_0
    }
    
    if(bestmse>innerfinalmse){
      bestmse=innerfinalmse
      bestf = override[nf]
      model.vault = list(newU=newU,beta_0=beta_0)
    }
    if(length(override)>1){if(nf==length(override)){print(paste('Final model has ',bestf,' factors.',sep=''))}}
    }
    grra$betas_grr = UtoB(newU=model.vault$newU,beta_0=model.vault$beta_0)$betas
    
  }
  lambda1_step2 =lambda1
  
  fcast <- c()
  if(TYPE==1){coef <- grr$betas_grr}
  if(TYPE==2){coef <- grra$betas_grr}

  #re-scale 
  for(j in 1:dim(X)[2]){
    coef[1,j+1,]=coef[1,j+1,]*scalingfactor[j]
  }
  coef[1,1,]=coef[1,1,]*sdy
  
  if(TYPE==1){grr$betas_grr<-coef
  grr$yhat =  sdy*grr$yhat}
  if(TYPE==2){grra$betas_grr<-coef
  grra$yhat =  sdy*grra$yhat}

  if(silent==0){
    plot.ts(cbind(msevec,diffvec))
    
    #what does beta look like
    par(mfrow=c(3,3))
    for(j in 1:dimX){ts.plot(betas_grr[1,j,])}
    par(mfrow=c(1,1))
    par(mfrow=c(3,3))
    for(j in 1:dimX){ts.plot(grra$betas_grr[1,j,])}
    par(mfrow=c(1,1))
  }
  
  if(length(oosX)>1){
    fcast <- t(as.matrix(append(1,oosX))) %*% as.matrix(coef[1,,dim(coef)[3]]) #basically the last beta_t
  }
  
  return(list(grr=grr,grrats=grrats,grra=grra,fcast=fcast,mse=mse,HPs=c(lambda1_step1,lambda1_step2,lambda2,r)))
}


gamma_reg_f = function(y,X,Z,adaptive=0,aparam=0,silent=1,lambdabooster=1,id=3,eweigths=1,alpha=0){
  Xplus = cbind(1,X)
  Zmod= Z
  #if(adaptive==1){Zmod=reweighter(Zmod,betas,aparam=aparam)}
  Zall = cbind(Xplus,Zmod)

  w=rep(1,dim(Zall)[2])
  w[1:dim(Xplus)[2]]=0 #1/dim(Zall)[1]

  #Transform y and x to impose identification
  r=dim(Zmod)[2]/dim(Xplus)[2]
  lambda_0 = matrix(0,dim(Xplus)[2],r)
  lambda_0[1:r,1:r]=diag(r)
  ytilde = y
  if(id==3){ytilde = y - Zmod%*%vec(t(lambda_0))}

  Zmod_idpart = Zmod
  if(id==3){
  selectmat = matrix(1,dim(Xplus)[2],r)
  selectmat[1:r,1:r]=0
  for(t in 1:dim(Zmod)[1]){Zmod_idpart[t,] = Zmod[t,]*vec(t(selectmat))}}
  #print((Zmod_idpart[1:4,]))
  Zall =  cbind(Xplus,Zmod_idpart)
  #ll = rep(-Inf,dim(Zall)[2])
  #ul = rep(Inf,dim(Zall)[2])
  #ll[(dim(Xplus)+1):(dim(Xplus)+r^2)]=vec(diag(r))
  #ul[(dim(Xplus)+1):(dim(Xplus)+r^2)]=vec(diag(r))
  
  #WILL NOT CHOSE THE PROPER LAMBDA FOR SAME REASONS AS ORIGINAL COSSO, should fix this at some point
  mse=NA
  if(lambdabooster!=0){CV_nnlasso = cv.glmnet(y=ytilde,x=Zall,intercept=0,alpha=alpha,
                         penalty.factor=w,standardize=FALSE,weights=eweigths,
                         lambda=exp(seq(log(0.001), log(5), length.out=100)))
  if(silent==0){plot(CV_nnlasso)}
  mse = min(CV_nnlasso$cvm)
  }
  
  lambdastar=0
  if(lambdabooster!=0){lambdastar=CV_nnlasso$lambda.min*lambdabooster}
  fit = glmnet(y=ytilde,x=Zall,lambda=lambdastar,standardize=FALSE,weights = eweigths,
              penalty.factor=w,intercept=0,alpha=alpha)
  gamma = as.vector(fit$beta)
  
  newLambda = t(matrix(gamma[(dim(Xplus)[2]+1):length(gamma)],nrow=r,ncol=dim(Xplus)[2]))
  
  if(id==3){newLambda[1:r,1:r]=diag(r)}
  
  return(list(mse=mse,gamma=gamma,Lambda=newLambda))
}

reweighterF = function(X,betas,aparam=0){
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