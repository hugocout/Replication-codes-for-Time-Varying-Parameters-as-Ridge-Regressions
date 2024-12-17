TVPRR_VARF <-function(X,Y,type,lambdavec=exp(linspace(-6,20,n=50)),lambda2=.1,
                 sweigths=1,oosX=1,kfold=5,silent=1,alpha=.5,nnlasso=0,olsprior=0,
                 adaptive=0,aparam=0,tol=10^(-2),maxit=15,starting_values=c(),CV.2SRR=TRUE,orthoFac=TRUE,
                 lambdabooster=1,var.share=.05,override=0,id=3,homo.param=.75,
                 sv.param=.75,fp.model=1,max.step.cv=1,CI=FALSE){
  TYPE=type
  #to do : 
  #-deal with a multivariate y
  #-do a loop of the first stage of TVP for each VAR eq (use parralel on this PC)
  #-factor analysis (putting all that togheter)
  #-stack Y in a vector, the rest should be standard but huge, 
    #if r smaller than M hopefully, solve the primal problem...actually if T*M > r*T + K*M 
  #-run un it on this PC
  
  ######### Data pre-processing ###################
  
  mse<-c() #
  grrats<-list()
  grra<- list()
  fcast<- c()
  
  #Standardize variances
  Y=as.matrix(Y)
  Yoriginal=as.matrix(Y)
  
  M=dim(Y)[2]
  bigT=dim(Y)[1]
  if(is.null(M)==TRUE){M=1}
  scalingfactor = array(1,dim=c(M,dim(X)[2]))
  sdy = 1:M
  
  for(j in 1:dim(X)[2]){
    for(m in 1:M){scalingfactor[m,j]=sd(Y[,m])/sd(X[,j])}
    X[,j]=(X[,j])/sd(X[,j])
  }
  for(m in 1:M){sdy[m]=sd(Y[,m])
  Y[,m] = Y[,m]/sdy[m]}

  ZZ <- Zfun(X) #Basis expansion
  YY <- as.matrix(Y) #univariate model
  yy = vec(YY)
  dimX= dim(as.matrix(X))[2]+1 # + constant
  
  if(is.null(starting_values)==TRUE){
  
  ######### TYPE 1 (FOR THE BASE MODEL OR ESTIMATION GIVEN HPS) ############
  #CV
    BETAS_GRR = array(0,dim=c(M,dimX,dim(Y)[1]))
    LAMBDAS = rep(0,M)
    EWmat = array(0,dim=dim(Y))
    YHAT_VAR=Y
    
    
    if(CV.2SRR==TRUE){
    if(length(lambdavec)>1){
      l.subset = 1:length(lambdavec) #c(round(length(lambdavec)/4),round(2*length(lambdavec)/4),round(3*length(lambdavec)/4))
      cvlist <- CV.KF.MV(ZZ=ZZ,Y=Y,dimX = dimX,lambdavec = lambdavec[l.subset],k=kfold,plot=abs(silent-1),sweigths=sweigths) #,eweigths = abs(yy)/mean(abs(yy)))
      lambdas_list =cvlist$lambdas_het
      }
    }
    else{
      lambdas_list = rep(M,lambdavec[round(length(lambdavec)/2)])
      mse=100
    }
    
  for(m in 1:M){
    #lambda1 <- cvlist$minimizer #lambdavec[round(length(lambdavec)/4)] 
    lambda1 <-  lambdas_list[m]
    mse <- 1000 #cvlist$minima
    
    #Final estimation
    grr<-dualGRR(ZZ,Y[,m],dimX,lambda1,olsprior = olsprior,lambda2=lambda2,sweigths=sweigths) #,eweigths = sqrt(abs(Y[,m]))/mean(abs(Y[,m])))
    
    BETAS_GRR[m,,] <- grr$betas_grr[1,,]
    LAMBDAS[m]=lambda1
    YHAT_VAR[,m]=grr$yhat
    
    e = Y[,m] - grr$yhat
    e[abs(e)>50*mad(e)]=0 #glicth detector
    #plot.ts(e)
    #garchstuff = hush(garchFit(data=e,~ garch(1,1)))
    #EW = (garchstuff@sigma.t)^sv.param
    #EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
    EWmat[,m] = 1 #EW
  }
  EWvec = vec(EWmat)/mean(EWmat)
  }
  else{ ####NEEDS TO BE UPDATED FOR MATRICES !!!!! #####
  BETAS_GRR=starting_values$BETAS_VARF_STD
  EWvec = starting_values$EWvec
  LAMBDAS=c()
  lambda1=starting_values$HPs[1]
  }
###### FACTOR TVP VAR #####
    grra <- list()
    maxit=maxit
    tol=tol

    msevec <- c() #mse1
    sigmavec <-c()
    diffvec <- c()
    new_lambdas_m=c()
    beta_0_w=rep(0,M)
    beta_0=vec(t(BETAS_GRR[,,1]))
    BETAS_GRRA = BETAS_GRR
    
    step=1
    par(mfrow=c(1,1))

    bestmse=Inf
    for(nf in 1:length(override)){
    innerbestmse=Inf
    
    if(sv.param>0){foolproof=2}
    else{foolproof=1}
    foolproof=1
    
    Xkron = kron(diag(M),cbind(1,X))
    
    BETAS_GRR_vec = matrix(BETAS_GRR,M*dimX,dim(Y)[1])
    aa = BtoF(BETAS_GRR_vec,var.share=var.share,override = override[nf],id=id,vizualizar = abs(silent-1),foolproof=foolproof)
    newF =aa$f
    newLambda = aa$Lambda
    
    if(is.null(starting_values)==FALSE){
    newF=starting_values$F
    newLambda=starting_values$Lambda
    }
      
    newU=newLambda%*%newF
    gamma = 1
    SW=1
    alpha_op=0
    
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
      NEWZ=c()
      for(m in 1:M){
        #newZ=R_Lambda2(X=X,Z=ZZ,newLambda = newLambda[((m-1)*dimX+1):(m*dimX),])
        newZ=R_Lambda(X=X,newLambda = newLambda[((m-1)*dimX+1):(m*dimX),])
        NEWZ=rbind(NEWZ,as.matrix(newZ))
      }
      newZkron = cbind(NEWZ[,1:(dim(newZ)[2]-dimX)],Xkron)
      yy_tilde = yy - Xkron%*%beta_0
      
      if(length(lambdavec)>1){if(-step>=-max.step.cv){
      lambda1pos = which(lambdavec==lambda1)
      l.subset = max(1,lambda1pos-4):min(length(lambdavec),lambda1pos+3)
      
      stuff <- CV.KF.MV(ZZ=newZkron,Y=yy_tilde,dimX = dimX*M,lambdavec = lambdavec[l.subset],
                                beta0_given = TRUE,k=kfold,plot=abs(silent-1),eweigths=EWvec,
                              sweigths=SW,nf=dim(newLambda)[2],lambda2=lambda2) #should it a be vec now

      lambda1=stuff$minimizer
      msevec = append(msevec,stuff$minima)
      }}
      else{
        lambda1=lambdavec
        msevec = append(msevec,0)
        }

      grra<-dualGRR(Zprime=newZkron,y=yy_tilde,dimX=dimX*M,lambda1=lambda1,eweigths=EWvec,
                    lambda2=lambda2,sweigths = SW,CI=0,olsprior=olsprior,
                    nf=dim(newLambda)[2])
      newF=grra$betas_grr[1,,2:dim(grra$betas_grr)[3]]-grra$betas_grr[1,,1:(dim(grra$betas_grr)[3]-1)]
      Yhatmat = matrix(grra$yhat,bigT,M)
      YBETA0mat=matrix(Xkron%*%beta_0,bigT,M)
      tot.Yhatmat = Yhatmat+YBETA0mat
      #ts.plot(Yhatmat)
      #ts.plot(YBETA0mat)
      #ts.plot(Y[,m])
      
      
      lambda2=1000000
      
      
      #SV update
      EWmat = array(1,dim=dim(Y))
      EWmat.std=EWmat
      if(sv.param>0){
      for(m in 1:M){
        e = Y[,m] - tot.Yhatmat[,m]
        e[abs(e)>10*mad(e)]=0 #glicth detector
        #plot(tot.Yhatmat[,m],Y[,m])
        
        garchstuff = hush(garchFit(data=e[2:length(e)],~ garch(1,1))) #arma(1,1) + 
        #svhat=exp(lm(log((e[2:length(e)]^2))~(log((e[1:(length(e)-1)]^2))))$fitted.values)
        #EW=append(svhat[1],svhat)
        EW = (garchstuff@sigma.t)^sv.param    #(svhat)^sv.param
        #EW = ((mean(garchstuff@sigma.t))^(1-sv.param))*(garchstuff@sigma.t)^sv.param    #(svhat)^sv.param
        
        EW=append(EW[1],EW)
        #if(-step>=-max.step.cv){
        #svstuff <- CV.KF.MV(ZZ=Zfun(as.matrix(tot.Yhatmat[,m]))[,c((bigT:(2*bigT-2)),2*bigT)],Y=yy_tilde,dimX = 1,lambdavec = lambdavec/100,
        #                  beta0_given = TRUE,k=kfold,plot=abs(silent-1),eweigths=1,
        #                  sweigths=1,nf=1,lambda2=0.00001) #should it a be vec now
        #}
        #lambda.sv =svstuff$minimizer
        #sv<-dualGRR(Zprime=Zfun(grra$yhat)[,c((bigT:(2*bigT-2)),2*bigT)],y=yy_tilde,dimX=1,lambda1=lambda.sv,eweigths=1,
        #             lambda2=0.00001,sweigths = 1,CI=0,olsprior=olsprior,
        #              nf=1)
        
        #EW=(hpfilter(e^2, freq=60,type="frequency")$trend)^sv.param
        #print(sv$betas_grr[1,,]
        #EW=(sv$betas_grr[1,,]) #^(2*sv.param)
        #print(EW)
        if(silent==0){plot.ts(EW)}
        
        #EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
        EWmat[,m] = EW
        EWmat.std[,m] = EW #/mean(EW)
        #EWmat.std[,m] = EWmat.std[,m]/EWmat.std[1,m]
        #EWmat[(length(e)-9):length(e),m]=exp(0.3*1:10)*EWmat[(length(e)-9):length(e),m]
      }}
      EWvec = vec(EWmat)/mean(EWmat)

      #Normalization
      utilT=dim(grra$betas_grr)[3]-1
      if(orthoFac==TRUE){
        #Normalization of F and creation of newZ for loadings regression
        if(length(newF)>utilT){
          orthomat = var(t(newF)) #(tcrossprod(newF)) #
          #SW = (diag(orthomat)/mean(diag(orthomat)))^(.5)
          orthoF= (pracma::sqrtm(orthomat)$Binv)%*%newF
        }
        else{orthoF=(newF)/sd(newF)} #-mean(newF)
        newF= orthoF
        #print(orthomat)
      }
      else{
        if(length(newF)>utilT){
          orthomat = diag(diag(var(t(newF)))) #(tcrossprod(newF)) #
          #SW = (diag(orthomat)/mean(diag(orthomat)))^(.5)
          orthoF= (pracma::sqrtm(orthomat)$Binv)%*%newF
        }
      }
      #else{for(jj in 1:r){newF[jj,]=(newF[jj,]-mean(newF[jj,]))/sd(newF[jj,])}} #}

      newZ= R_f2(f=newF,dimX=dimX,Z=ZZ)
      
      newLambda=c()
      beta_0=c()
      lambda.glmnet=c()
      old_lambdas_m=new_lambdas_m
      if(M>4){par(mfrow=c(3,3))}
      else{par(mfrow=c(2,2))}
      
      if(step>0.5*max.step.cv){alpha_op=alpha}
      
      
      for(m in 1:M){
      if(step>max.step.cv){lambda.glmnet=old_lambdas_m[m]}
      gammareg= gamma_reg_f_VAR(y=Y[,m],X=cbind(X),Z=newZ,eweigths=EWmat.std[,m],alpha=alpha_op,
                                silent=1,lambdabooster=lambdabooster*min(1,10*((step/(max.step.cv-1))^2)),impose.lambda=lambda.glmnet,beta_0_w=beta_0_w[m])
      gamma_new=gammareg$gamma
      new_lambdas_m=append(new_lambdas_m,gammareg$lambdastar)
      beta_0 = append(beta_0,gamma_new[1:(dimX)])
      if(-step>=-2*max.step.cv){beta_0_w[m] = 0.01+mean(abs(gammareg$Lambda)+0.01)/mean(abs(gamma_new[1:(dimX)])+0.01)}
      gamma_new = gamma_new[(dimX+1):(length(gamma_new))] 
      gamma=gamma_new
      if(silent==0){barplot(t(gamma))}
      newLambda = rbind(newLambda,gammareg$Lambda)
      }
      #print(mean(abs(gammareg$Lambda)))
      #print(mean(abs(beta_0)))

      par(mfrow=c(1,1))
    
      if(r>foolproof){
      dumpvec=c()
      for(rr in 1:r){if(-sum(abs(newLambda[,rr]))>-0.0001*tol){dumpvec=append(dumpvec,rr)}}
      if(is.null(dumpvec)==FALSE){
        newLambda=as.matrix(newLambda[,-dumpvec])
        newF=as.matrix(newF[-dumpvec,])
        if(dim(newF)[1]>dim(newF)[2]){newF=t(newF)}
        r = r - length(dumpvec)
      }}
      else{if(-sum(abs(newLambda))>-0.01*tol){r=0}}
      
      #Impose normalization of Lambda matrix
      oldU=newU
      newU=newLambda%*%newF
      
      if(r==0){
        print(paste('converged in', step,'steps because all variances are 0'))
        break
      }
      
      if(id==1){
      #if(step > 2*max.step.cv){var.share=0}
        print(var.share*min(1,sqrt(2*step/max.step.cv)))
      aa = UtoF((newU),var.share=var.share*min(1,sqrt(2*step/max.step.cv)),override = r,id=id,vizualizar = abs(silent-1),foolproof=foolproof)
      if(orthoFac==TRUE){
      newF =aa$f
      newLambda = aa$Lambda
      }
      if(r!=dim(newLambda)[2]){
        maxit =maxit+5
        SW=1}
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
    #print(beta_0)
    coef = UtoB(newU=model.vault$newU,beta_0=model.vault$beta_0)$betas

    f_CI=c()
    if(CI==TRUE){
      
      f_CI<-dualGRR(Zprime=newZkron,y=yy_tilde,dimX=dimX*M,lambda1=lambda1,eweigths=EWvec,
                  lambda2=lambda2,sweigths = SW,CI=1.96,olsprior=olsprior,
                  nf=dim(newLambda)[2])$betas_grr_ci
      
      #for(m in 1:M){
      #  lambda.glmnet=old_lambdas_m[m]
      #  gammareg= gamma_reg_f_VAR(y=Y[,m],X=cbind(1,X),Z=newZ,eweigths=EWmat.std[,m],alpha=alpha_op,
      #                            silent=1,lambdabooster=lambdabooster,impose.lambda=lambda.glmnet,beta_0_w=beta_0_w[m])
      #  gamma_new=gammareg$gamma
      #  new_lambdas_m=append(new_lambdas_m,gammareg$lambdastar)
      #  beta_0_w[m] 
      #  newLambda_ci[,,1] = rbind(newLambda,gammareg$Lambda)
      #  newLambda_ci[,,2] = rbind(newLambda,gammareg$Lambda)
      #}
    #}
      
      
      
      }
    
    
  lambda1_step2 =lambda1
  
  fcast <- c()
  #Create matrix version
  BETAS_VARF=BETAS_GRR
  for(m in 1:M){BETAS_VARF[m,,]=coef[1,((m-1)*dimX+1):(m*dimX),]}
  BETAS_VARF_STD=BETAS_VARF
  RESIDUALS_STD=Y-tot.Yhatmat
  if(is.null(starting_values)==TRUE){RESIDUALS_VAR_STD=Y-YHAT_VAR}
  else{RESIDUALS_VAR_STD=c()}
  
  
  #re-scale 
  for(m in 1:M){for(j in 1:(dimX-1)){
    BETAS_VARF[m,j+1,]=BETAS_VARF[m,j+1,]*scalingfactor[m,j]
    BETAS_GRR[m,j+1,]=BETAS_GRR[m,j+1,]*scalingfactor[m,j]
  }
    BETAS_VARF[m,1,]=BETAS_VARF[m,1,]*sdy[m]
    BETAS_GRR[m,1,]=BETAS_GRR[m,1,]*sdy[m]
    tot.Yhatmat[,m]=tot.Yhatmat[,m]*sdy[m]
  }
  RESIDUALS=Yoriginal-tot.Yhatmat
  if(is.null(starting_values)==TRUE){RESIDUALS_VAR=Yoriginal-YHAT_VAR}
  else{RESIDUALS_VAR=c()}
  
  #coef[1,1,]=coef[1,1,]*sdy
  #coef[1,j+1,]=coef[1,j+1,]*scalingfactor[m,j]

  #if(silent==0){
  #  plot.ts(cbind(msevec,diffvec))
  #  #what does beta look like
  #  par(mfrow=c(3,3))
  #  for(m in 1:M){for(j in 1:(dimX)){ts.plot(BETAS_VARF[m,j,])}}
  #  par(mfrow=c(1,1))
  #}
  
  if(length(oosX)>1){
    #OOSX = o.call(rbind, replicate(M, as.matrix(append(1,oosX)), simplify=FALSE))
    fcast <-  BETAS_VARF[,,dim(BETAS_VARF)[3]] %*% as.matrix(append(1,oosX))  #basically the last beta_t
  }
  
  HPs=c(lambda1,lambda2,r)
  return(list(BETAS_VAR=BETAS_GRR,BETAS_VARF=BETAS_VARF,YHAT=tot.Yhatmat,f_CI=f_CI,
              LAMBDAS=LAMBDAS,grra=grra,fcast=fcast,mse=mse,HPs=HPs,
              res=list(RESIDUALS_VARF=RESIDUALS,RESIDUALS_VARF_STD=RESIDUALS_STD,
                       RESIDUALS_VAR=RESIDUALS_VAR,RESIDUALS_VAR_STD=RESIDUALS_VAR_STD),
              starter_pack=list(EWvec=EWvec,F=newF,Lambda=newLambda,BETAS_VARF_STD=BETAS_VARF_STD,HPs=HPs)))
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