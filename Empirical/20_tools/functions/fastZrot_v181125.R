BtoF <- function(betas,var.share=.8,override=0,vizualizar=0,id=1,foolproof=2){
  
  lbeta=dim(betas)[2]
  t=lbeta
  
  umat = betas[,2:lbeta]-betas[,1:(lbeta-1)]
  omega = umat%*%t(umat)
  ev<- eigen(omega)$values/sum(eigen(omega)$values)
  normev = ev/max(ev)
  cev = cumsum(ev)
  if(vizualizar==1){barplot(normev)}
  
  p= dim(umat)[1]
  fanal = factor(t(umat),dim(umat)[1]) #factor analysis 
  #print(mean(diag(omega)))
  
  #print(mean(diag(t(fanal$factors)%*%fanal$factors)))
  newev = array(0,dim=c(1,p)) 
  evnotzero = eigen(omega)$values[normev>var.share] #(selection rule for now) no seledtion here
  #evnotzero = eigen(omega)$values[(cev)<var.share] #(selection rule for now) no seledtion here
  if(length(evnotzero)>0){newev[1,1:length(evnotzero)] =  evnotzero} #put zeros on the other eigenvalues
  else{newev[1,1] =  eigen(omega)$values[1]}
  nf = length(evnotzero)
  
  #Bai ng tentative implementation
  #################
  #baing=0
  #if(baing==1){
  #ICvec=c()
  #for(j in 1:round(p-1)){
  #  ICvec[j]=ICp2(mse=factor(t(umat),n_fac=j,id=1)$mse,n_fac=j, bigt=t, bign=p)
  #}
  #if(vizualizar==1){barplot(ICvec)}
  #nf = which(ICvec==min(ICvec))
  #}
  nf=max(min(override,nf),foolproof) #the last max is the foolproof for ID
  ##################
    
  print(paste(nf,'factors selected'))
  
  fanal = factor(t(umat),nf,id=id) #factor analysis
  
  return(list(f=t(as.matrix(fanal$factors)),Lambda=fanal$lambda))
}


UtoF <- function(umat,var.share=.8,override=0,vizualizar=0,id=1,foolproof=2){
  
  omega = umat%*%t(umat)
  ev<- eigen(omega)$values/sum(eigen(omega)$values)
  normev = ev/max(ev)
  cev = cumsum(ev)
  if(vizualizar==1){barplot(normev)}
  
  p= dim(umat)[1]
  fanal = factor(t(umat),dim(umat)[1]) #factor analysis 
  #print(mean(diag(omega)))
  
  #print(mean(diag(t(fanal$factors)%*%fanal$factors)))
  newev = array(0,dim=c(1,p)) 
  evnotzero = eigen(omega)$values[normev>var.share] #(selection rule for now) no seledtion here
  #evnotzero = eigen(omega)$values[(cev)<var.share] #(selection rule for now) no seledtion here
  if(length(evnotzero)>0){newev[1,1:length(evnotzero)] =  evnotzero} #put zeros on the other eigenvalues
  else{newev[1,1] =  eigen(omega)$values[1]}
  nf =length(evnotzero) #ID foolproof
  
  #################
  #baing=0
  #if(baing==1){
  #  ICvec=c()
  #  for(j in 1:round(p/1)){
  #    ICvec[j]=ICp2(mse=factor(t(umat),n_fac=j,id=1)$mse,n_fac=j, bigt=t, bign=p)
  #  }
  #  if(vizualizar==1){barplot(ICvec)}
  #  nf = which(ICvec==min(ICvec))
  #}
  nf=max(min(override,nf),foolproof)
  ################
  
  #print(paste(nf,'factors selected'))
  
  fanal = factor(t(umat),nf,id=id) #factor analysis
  
  return(list(f=t(as.matrix(fanal$factors)),Lambda=fanal$lambda))
}

R_Lambda3<- function(X,newLambda){
  
  r=dim(newLambda)[2]
  if(is.null(r)==TRUE){r=1}
  
  XX= cbind(matrix(1,(nrow(X)),1),X)
  Z=array(0,dim=c((nrow(X)),(nrow(X)-1),r))
  
  #The sequence of tau's
  param_seq  <- seq(2,nrow(X),by=1)
  
 # if(r!=1){
  #New regressors
  for(tt in 2:(dim(X)[1])){
        for(jj in 1:r){
          Z[tt,1:(tt-1),jj]= repmat(as.matrix(newLambda)[,jj]%*%XX[tt,],tt-1,1)
        }
 
  }
  Zprime <- c()
  for(kk in 1:r){
      Zprime = cbind(Zprime,Z[,,kk])
    }  

  Zrot <- cbind(Zprime,XX)
  return(Zrot)
}

R_Lambda<- function(X,newLambda){
  
  r=dim(newLambda)[2]
  if(is.null(r)==TRUE){r=1}
  
  XX= cbind(matrix(1,(nrow(X)),1),X)
  Z=array(0,dim=c((nrow(X)),(nrow(X)-1),r))
  
  #New regressors
  for(tt in 2:(dim(X)[1])){
    #for(tho in 2:tt){
    #if(tho<=tt){
    Z[tt,1:(tt-1),]= repmat(t(crossprod(as.matrix(newLambda)[,],XX[tt,])),tt-1,1)
    
  }
  Zprime <- c()
  for(kk in 1:r){
    Zprime = cbind(Zprime,Z[,,kk])
  }  
  
  
  Zrot <- cbind(Zprime,XX)
  return(Zrot)
}

R_Lambda2<- function(Z,X,newLambda){
  
  r=dim(newLambda)[2]
  if(is.null(r)==TRUE){r=1}
  
  Xnew= cbind(matrix(1,(nrow(X)),1),X)
  dimX=dim(Xnew)[2]
  
  #print(dim(Z[,1:(dim(Z)[2]-dimX)]))
  #print(dim(kron(newLambda,diag(dim(Xnew)[1]-1))))
  
  Zprime <-Z[,1:(dim(Z)[2]-dimX)]%*% kron(newLambda,diag(dim(Xnew)[1]-1))
  
  Zrot <- cbind(Zprime,Xnew)
  return(Zrot)
}

R_f<- function(f,X){
  
  f = t(as.matrix(f))
  r= dim(f)[2]
  
  #f = f - mean(f) #identification (factors must have mean 0)
  
  #Good old R vector matrix problem
  if(r==nrow(X)-1){
      f=t(f)
      r= dim(f)[2]
  }
  
  nsF=f
  #nsF[1,]=f[1,]
  
  #canceled this, have to normalize in a multivariate way
  #for(j in 1:r){
  #  f[,j]=(f[,j])/sd(f[,j])
  #}
  
  #print((dim(f)))
  #for(tt in 2:(dim(f)[1])){      # MODIFY DUALLGRR ACCORDINGLY!!!
  #    nsF[tt,] = f[tt,]+nsF[tt-1,]
  #}
  #new = lm(nsF ~ seq(1,length(nsF)))$residuals
  #std = 0.8/(range(new)[2]-range(new)[1])
  #nsF = as.matrix(std*new)
  
  Zrot = c()
  Xnew= cbind(matrix(1,(nrow(X)),1),X)

  for(j in 1:r){
    for(k in 1:dim(Xnew)[2]){
      Xf = Xnew[,k]*nsF[,j]
      Zrot = cbind(Zrot,Xf)
    }
  }
  #ts.plot(XX[,1])
  #ts.plot(nsF[,1])
  #W=array(0,dim=c((nrow(X)),(nrow(X)-1),(ncol(X)+1)))
  
  #The sequence of tau's
  #param_seq  <- seq(2,nrow(X),by=1)
  
  #New regressors
  #  for(tho in param_seq){W[tho,tho-1,]= XX[tt,]}
  
  #Wprime <- c()
  #for(kk in 1:(ncol(X)+1)){
  #  Wprime = cbind(Wprime,W[,,kk])
  #}
  
  #Zprime <- Wprime[,]%*%(kron(nsF,diag(dim(X)[2]+1)))   #[,1:r] t()

  #Zrot <- Zprime #cbind(Zprime,XX)
  return(Zrot)
}

R_f2 = function(f,Z,dimX){
  
  #Good old R vector matrix problem
  r= dim(f)[2]
  if(is.null(r)==TRUE){f=t(as.matrix(f))}
  
  Zrot = Z[,1:(dim(Z)[2]-dimX)] %*% kron(diag(dimX),t(f))
  return(Zrot)
}

FtoB <- function(Lambda,f,beta_0){
  
  r = dim(f)[2]
  dimX=length(beta_0)
  F = f[1,,2:dim(f)[3]]-f[1,,1:(dim(f)[3]-1)]

  #Make U matrix
  if(r!=1){newU = Lambda%*%F}
  else{newU = (Lambda)%*%t(F)}

  #Generate final betas
  eq=1
  newbetas = array(0,dim=c(1,dimX,dim(f)[3]))
  newbetas[1,,1]=beta_0
  
  for(tt in 2:(dim(f)[3])){
    newbetas[eq,,tt] = newbetas[eq,,tt-1] + newU[,(tt-1)]
  }
  
  return(list(betas=newbetas,nsF=f))
}

UtoB <- function(newU,beta_0){
  
  #Generate final betas
  eq=1
  newbetas = array(0,dim=c(1,length(beta_0),dim(newU)[2]+1))
  newbetas[1,,1]=beta_0
  #+newU[,1]
  
  
  for(tt in 2:(dim(newU)[2]+1)){
    newbetas[eq,,tt] = newbetas[eq,,tt-1] + newU[,(tt-1)]
  }
  
  newbetas[1,,1]=beta_0+newU[,1]
  
  return(list(betas=newbetas))
}

FtoB2 <- function(Lambda,f,beta_0){
  
  r = dim(f)[2]
  dimX=length(beta_0)
  F = f[1,,2:dim(f)[3]] #-f[1,,1:(dim(f)[3]-1)]
  
  #Make U matrix
  if(r!=1){newU = Lambda%*%F}
  else{newU = (Lambda)%*%t(F)}
  
  #Generate final betas
  eq=1
  newbetas = array(0,dim=c(1,dimX,dim(f)[3]))
  newbetas[1,,1]=beta_0
  
  for(tt in 2:(dim(f)[3])){
    newbetas[eq,,tt] = newU[,(tt-1)] # newbetas[eq,,tt-1] +
  }
  
  return(list(betas=newbetas,nsF=f))
}

ICp2 <- function(mse, n_fac, bigt, bign){
  
  k <- n_fac
  v <- mse
  
  c_nt <- min(c(sqrt(bign), sqrt(bigt)))
  CT <- k*((bign+bigt)/(bign*bigt))*log(c_nt^2) # Penalty function
  
  icp2 <- log(v) + CT
  
  return(icp2)
}

gamma_reg_f_VAR = function(y,X,Z,silent=1,lambdabooster=1,
                           eweigths=1,alpha=0,kfold=5,impose.lambda=c(),beta_0_w=0){
  Xplus = X
  dimXplus=dim(Xplus)[2]+1
  Zmod= Z
  Zall = cbind(Xplus,Zmod)
  
  w=rep(1,dim(Zall)[2])
  w[1:dimXplus]=min(100,max(0.01,beta_0_w)) #1/dim(Zall)[1]
  
  #Transform y and x to impose identification
  r=dim(Zmod)[2]/dimXplus
  
  #WILL NOT CHOSE THE PROPER LAMBDA FOR SAME REASONS AS ORIGINAL COSSO, should fix this at some point
  mse=NA
  lambdastar=0
  if(lambdabooster!=0){if(is.null(impose.lambda)==TRUE){
    
    CV_nnlasso = cv.glmnet(y=y,x=Zall,intercept=TRUE,alpha=alpha,
                           penalty.factor=w,standardize=FALSE,weights=eweigths,nfolds=5,
                           lambda=exp(seq(log(0.001), log(9), length.out=15)))
    if(silent==0){plot(CV_nnlasso)}
    mse = min(CV_nnlasso$cvm)
    lambdastar=CV_nnlasso$lambda.min*lambdabooster
  }
    else{lambdastar=impose.lambda}
  }
  
  fit = glmnet(y=y,x=Zall,lambda=lambdastar,standardize=FALSE,weights = eweigths,
               penalty.factor=w,intercept=TRUE,alpha=alpha)
  
  gamma = append(fit$a0,as.vector(fit$beta))
  
  newLambda = t(matrix(gamma[(dimXplus+1):length(gamma)],nrow=r,ncol=dimXplus))
  #print(gamma)
  return(list(mse=mse,gamma=gamma,Lambda=newLambda,lambdastar=lambdastar))
}