# =============================================================================================================
# INITIALIZATION
# =============================================================================================================

rm(list=ls()) # To clear the current environment
yymmdd = paste(substr(Sys.Date(),3,4),substr(Sys.Date(),6,7),substr(Sys.Date(),9,10),sep='')

wd = 'INSERT YOUR PATH HERE/Empirical/'
setwd(wd)
paths <- list(pro = "00_prog",
              dat = "10_data",
              too = "20_tools",
              fun = "20_tools/functions",
              out = "30_output",
              rst = "40_results")

#Install (for Cluster only)
# Array of packages to be checked and installed if not already
myPKGs <- c('readr', 'GA', 'e1071', 'pracma','glmnet','fGarch','matrixcalc')

# check to see if each package is installed, and if not add it to a list to install
InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]

# install any needed packages
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

# Libraries
library(readr)      # To import data
library(GA)         # For genetic algorithm
library(e1071)      # To estimate SVR
library(pracma)     # For repmat 
library(glmnet)
library(fGarch)
library(matrixcalc)

# Custom Functions -----------------------------------------------
source(paste(paths$too, 'EM_sw.R', sep='/'))
source(paste(paths$too, 'ICp2.R', sep='/'))
source(paste(paths$too, 'Xgenerators_v190127.R', sep='/'))

# Functions for TVPRR
source(paste(paths$fun, 'dualGRRmdA_v190215.R', sep='/'))
source(paste(paths$fun, 'CVGSBHK_v181127.R', sep='/'))
source(paste(paths$fun, 'CVKFMV_v190214.R', sep='/'))
source(paste(paths$fun, 'MV2SRR_v200722.R', sep='/'))

source(paste(paths$fun, 'zfun_v180722.R', sep='/'))
source(paste(paths$fun, 'factor.R', sep='/'))
source(paste(paths$fun, 'TVPRRcosso_v181120.R', sep='/'))
source(paste(paths$fun, 'TVPRRcossoF_v190125.R', sep='/'))
source(paste(paths$fun, 'TVPRR_v181111.R', sep='/'))
source(paste(paths$fun, 'fastZrot_v181125b.R', sep='/'))
source(paste(paths$fun, 'TVPRR_VARF_v190304.R', sep='/'))

# =============================================================================================================
# DATA IMPORT
# =============================================================================================================

dates_all          <- seq(as.Date('1973-01-01'), 
                          as.Date('2016-01-01'), 
                          by = "month")
cs18_f4         <- as.matrix(read_csv(paste(paths$dat,'cs18_fig4_ppuzzle.csv',sep='/')))

# =============================================================================================================
# ESTIMATION
# =============================================================================================================

gap=1
lambdavec=exp(linspace(4,18,n=15))
silenceplz = 1
scree_ts=.02

for(v in 1:3){
  cBETAS=array(NA,dim=c(48,2,3,419))
  mod = 2
  
  ################# GDP ################
  if(v==1){
    depvars = 97:(97+47)
    Ymat=cs18_f4[,depvars]
    X = cs18_f4[,1:96]  #65 is 48 lags of shock...let's begin with 12...30
    szv = dim(Ymat)[2]
    train = cbind(Ymat,X)
    
    minnna = rep(NA,dim(train)[2])
    maxnna=minnna
    for(m in 1:dim(train)[2]){
      minnna[m] <- min(which(!is.na(train[,m])))
      maxnna[m] <- max(which(!is.na(train[,m])))
    }
    subset = max(minnna):min(maxnna)
    date_fig4=dates_all[subset]
    lp.b1.pos=7
  }
  
  ################# UR ################
  if(v==2){
    depvars=c(97,97+48,97+48*2) 
    depvars = (97+48+48):(97+48+48+47)
    Ymat=cs18_f4[,depvars]
    X = cs18_f4[,1:96]  #65 is 48 lags of shock...let's begin with 12...30
    szv = dim(Ymat)[2]
    train = cbind(Ymat,X)
    
    minnna = rep(NA,dim(train)[2])
    maxnna=minnna
    for(m in 1:dim(train)[2]){
      minnna[m] <- min(which(!is.na(train[,m])))
      maxnna[m] <- max(which(!is.na(train[,m])))
    }
    subset = max(minnna):min(maxnna)
    date_fig4=dates_all[subset]
    lp.b1.pos=19
    
  }
  
  ################# INF ################
  if(v==3){
    depvars=c(97,97+48,97+48*2) 
    depvars = (97+48):(97+48+47)
    Ymat=cs18_f4[,depvars]
    X = cs18_f4[,1:96] 
    szv = dim(Ymat)[2]
    train = cbind(Ymat,X)
    # 1. Plain
    
    minnna = rep(NA,dim(train)[2])
    maxnna=minnna
    for(m in 1:dim(train)[2]){
      minnna[m] <- min(which(!is.na(train[,m])))
      maxnna[m] <- max(which(!is.na(train[,m])))
    }
    subset = max(minnna):min(maxnna)
    date_fig4=dates_all[subset]
    lp.b1.pos=1
  }
  
  ##############################################
  ############# ESTIMATION #####################
  ##############################################
  cvvec = c()
  for(h in 1:48){ 
    CV = cv.glmnet(x=X[subset,],y=Ymat[subset,h],  family='gaussian', alpha=0)
    cvvec=append(cvvec,CV$lambda.min)
  }
  
  if(mod==1){ 
    gdp1 <- MV.2SRR(Y=Ymat[subset,],X=X[subset,],orthoFac = TRUE,CV.2SRR = TRUE,CV.again = FALSE,
                    lambdavec =lambdavec,sweigths=1,type=2,fp.model=1,sv.param = 0,homo.param = 0,
                    alpha=0,silent=0*silenceplz,kfold=5,lambda2=mean(cvvec)/2,max.step.cv = 5,
                    adaptive=1,aparam=-.5,tol=10^(-10),maxit=30,starting_values=c(), #lp.b1.pos=lp.b1.pos,
                    lambdabooster=1,var.share = scree_ts,override = 5,id=1,oosX=c(),CI=TRUE)
  }
  
  if(mod==2){ #-c(25:48)
    gdp1 <- MV.2SRR(Y=Ymat[subset,],X=X[subset,],orthoFac = TRUE,CV.2SRR = TRUE,CV.again =TRUE,
                    lambdavec =lambdavec,sweigths=1,type=2,fp.model=1,sv.param = .75,homo.param = .75,
                    alpha=0,silent=0*silenceplz,kfold=5,lambda2=mean(cvvec)/2,max.step.cv = 5,
                    adaptive=1,aparam=-.5,tol=10^(-10),maxit=30,starting_values=c(), #lp.b1.pos=lp.b1.pos,
                    lambdabooster=1,var.share = scree_ts,override = 5,id=1,oosX=c(),CI=TRUE)
  }
  
  
  for(t in 1:419){
    cBETAS[,mod,v,t]=cumsum(gdp1$BETAS_2SRR[,50,t]) 
  }
  
  BOLS=solve(crossprod(cbind(1,X[subset,])))%*%crossprod(cbind(1,X[subset,]),Ymat[subset,])
  YHAT_OLS = cbind(1,X[subset,])%*%BOLS #-c(25:48)
  
  filename= paste(paths$out,'/figure3/cs18results',v,mod,'.Rdata',sep='')
  save(gdp1,cBETAS,BOLS,YHAT_OLS,mod,v,file=filename)
}
