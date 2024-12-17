MV.2SRR <-function(X,Y,type,lambdavec=exp(linspace(-6,20,n=50)),lambda2=.1,CV.again=FALSE,
                      sweigths=1,oosX=1,kfold=5,silent=1,alpha=.5,nnlasso=0,olsprior=0,
                      adaptive=0,aparam=0,tol=10^(-2),maxit=15,starting_values=c(),CV.2SRR=TRUE,orthoFac=TRUE,
                      lambdabooster=1,var.share=.05,override=0,id=3,homo.param=.75,
                      sv.param=.75,fp.model=1,max.step.cv=1,CI=FALSE){
  TYPE=type

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
  
 # if(is.null(starting_values)==TRUE){
    
    ######### TYPE 1 (FOR THE BASE MODEL OR ESTIMATION GIVEN HPS) ############
    #CV
    BETAS_GRR = array(0,dim=c(M,dimX,dim(Y)[1]))
    BETAS_GRRATS =BETAS_GRR
    LAMBDAS = rep(0,M)
    EWmat = array(0,dim=dim(Y))
    YHAT_VAR=Y
    YHAT_VARF=Y
    
    
    
    if(CV.2SRR==TRUE){
      if(length(lambdavec)>1){
        l.subset = 1:length(lambdavec) #c(round(length(lambdavec)/4),round(2*length(lambdavec)/4),round(3*length(lambdavec)/4))
        cvlist <- CV.KF.MV(ZZ=ZZ,Y=Y,dimX = dimX,lambdavec = lambdavec[l.subset],lambda2=lambda2,
                           k=kfold,plot=abs(silent-1),sweigths=sweigths) #,eweigths = abs(yy)/mean(abs(yy)))
        lambdas_list =cvlist$lambdas_het
      }
    }
    else{
      lambdas_list = rep(M,lambdavec[round(length(lambdavec)/2)])
      mse=100
    }
    
    for(m in 1:M){
      print(m)
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
      garchstuff = hush(garchFit(data=e,~ garch(1,1)))
      EW = (garchstuff@sigma.t)^sv.param
      EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
      EWmat[,m] = EW
      
      #NEW ###############
      betas_grr = BETAS_GRR[m,,]
      umat = betas_grr[,2:nrow(ZZ)]-betas_grr[,1:(nrow(ZZ)-1)]
      sigmasq= diag(umat%*%t(umat))^homo.param
      sigmasq=(sigmasq/mean(sigmasq))
      
      #e = yy - grr$yhat
      #EW = (hush(garchFit(data=e,~ garch(1,1)))@sigma.t)^sv.param
      #EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement 
      
      if(CV.again){cvlist <- cvgs.bhk2015(ZZ=ZZ,Y=Y[,m],dimX = dimX,lambda2 = lambda2,eweigths = EW,
                                                     lambdavec = lambdavec,k=kfold,plot=abs(silent-1),sweigths=sigmasq)
      usethislambda1=cvlist$minimizer
      }
      else{usethislambda1=lambda1}
      
      if(homo.param>0 | sv.param>0){
      grrats<-dualGRR(ZZ,Y[,m],dimX,lambda1=usethislambda1,eweigths = EW,olsprior = olsprior,
                      lambda2=lambda2,sweigths = sigmasq)
      betas_grrats <- grrats$betas_grr
      BETAS_GRRATS[m,,]=grrats$betas_grr[1,,]
      YHAT_VARF[,m]=grrats$yhat
      
      }else{
        betas_grrats <- grr$betas_grr
        BETAS_GRRATS[m,,]=grr$betas_grr[1,,]
        YHAT_VARF[,m]=YHAT_VAR[,m]
      }
      
      
    }
    EWvec = vec(EWmat)/mean(EWmat)
  #}
  # else{ ####NEEDS TO BE UPDATED FOR MATRICES !!!!! #####
  #   BETAS_GRR=starting_values$BETAS_VARF_STD
  #   EWvec = starting_values$EWvec
  #   LAMBDAS=c()
  #   lambda1=starting_values$HPs[1]
  # }

  
  lambda1_step2 =lambda1
  
  fcast <- c()
  #Create matrix version
  BETAS_VARF=BETAS_GRRATS
  #for(m in 1:M){BETAS_VARF[m,,]=coef[1,((m-1)*dimX+1):(m*dimX),]}
  BETAS_VARF_STD=BETAS_VARF
  #RESIDUALS_STD=Y-tot.Yhatmat
  #if(is.null(starting_values)==TRUE){RESIDUALS_VAR_STD=Y-YHAT_VAR}
  #else{RESIDUALS_VAR_STD=c()}
  
  
  #re-scale 
  for(m in 1:M){for(j in 1:(dimX-1)){
    BETAS_VARF[m,j+1,]=BETAS_VARF[m,j+1,]*scalingfactor[m,j]
    BETAS_GRR[m,j+1,]=BETAS_GRR[m,j+1,]*scalingfactor[m,j]
  }
    BETAS_VARF[m,1,]=BETAS_VARF[m,1,]*sdy[m]
    BETAS_GRR[m,1,]=BETAS_GRR[m,1,]*sdy[m]
    YHAT_VARF[,m]= YHAT_VARF[,m]*sdy[m]
    YHAT_VAR[,m]= YHAT_VAR[,m]*sdy[m]
    
  }
  #RESIDUALS=Yoriginal-tot.Yhatmat
 # if(is.null(starting_values)==TRUE){RESIDUALS_VAR=Yoriginal-YHAT_VAR}
  #else{RESIDUALS_VAR=c()}

  
  if(length(oosX)>1){
    #OOSX = o.call(rbind, replicate(M, as.matrix(append(1,oosX)), simplify=FALSE))
    fcast <-  BETAS_VARF[,,dim(BETAS_VARF)[3]] %*% as.matrix(append(1,oosX))  #basically the last beta_t
  }
  
  HPs=c(lambda1,lambda2)
  return(list(BETAS_RR=BETAS_GRR,BETAS_2SRR=BETAS_VARF,#f_CI=f_CI, #YHAT=tot.Yhatmat,
              LAMBDAS=LAMBDAS,grra=grra,fcast=fcast,mse=mse,HPs=HPs,
              YHAT_RR=YHAT_VAR,YHAT_2SRR=YHAT_VARF,Y=Yoriginal,
              #res=list(RESIDUALS_VARF=RESIDUALS,RESIDUALS_VARF_STD=RESIDUALS_STD,
               #        RESIDUALS_VAR=RESIDUALS_VAR,RESIDUALS_VAR_STD=RESIDUALS_VAR_STD),
              starter_pack=list(EWvec=EWvec,BETAS_VARF_STD=BETAS_VARF_STD,HPs=HPs)))
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