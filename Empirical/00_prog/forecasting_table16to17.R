
# =============================================================================================================
# INITIALIZATION
# =============================================================================================================

rm(list=ls()) # To clear the current environment
wd = 'INSERT YOUR PATH HERE/Empirical/'
setwd(wd)

paths <- list(pro = "00_prog",
              dat = "10_data",
              too = "20_tools",
              fun = "20_tools/functions",
              out = "30_output",
              rst = "40_results")

# Array of packages to be checked and installed if not already
myPKGs <- c('readr', 'GA', 'e1071', 'pracma','doParallel','foreach','glmnet','timeSeries','fGarch','matrixcalc')

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
library(doParallel) # For detectCores
library(foreach)    # For foreach
library(glmnet)
library(timeSeries)
library(fGarch)
library(matrixcalc)


# =============================================================================================================
# FORECASTING PARAMETERS
# =============================================================================================================

# Create all possible combinations for model selection, variables and horizons ------------------------------
combn <- list(M = c(1:4),
              V = c(1:5),
              H = c(1,2,4))
all_options <- expand.grid(combn)
rownames(all_options) <- c()

# Number of cores to use (all but one, if set to 1 there is no parallelization)
ncores <- detectCores() - 1

# Forecasting function ---------------------------------------------------------------------------------------
POOS <- function(it_pos) {
  
  # Libraries
  library(readr)      # To import data
  library(GA)         # For genetic algorithm
  library(e1071)      # To estimate SVR
  library(pracma)     # For repmat 
  library(doParallel) # For detectCores
  library(foreach)    # For foreach
  library(glmnet)
  library(timeSeries)
  library(fGarch)
  library(matrixcalc)
  
  # General
  source(paste(paths$too, 'EM_sw.R', sep='/'))
  source(paste(paths$too, 'ICp2.R', sep='/'))
  source(paste(paths$too, 'Xgenerators_v190127.R', sep='/'))
  
  # Functions for TVPRR
  source(paste(paths$fun, 'dualGRRmdA_v190215.R', sep='/'))
  source(paste(paths$fun, 'CVGSBHK_v181127.R', sep='/'))
  source(paste(paths$fun, 'zfun_v190304.R', sep='/'))
  source(paste(paths$fun, 'factor.R', sep='/'))
  source(paste(paths$fun, 'TVPRRcosso_v181120.R', sep='/'))
  source(paste(paths$fun, 'TVPRRcossoF_v190125.R', sep='/'))
  source(paste(paths$fun, 'TVPRR_v181111.R', sep='/'))
  source(paste(paths$fun, 'fastZrot_v181125.R', sep='/'))
  source(paste(paths$fun, 'CVKFMV_v190214.R', sep='/'))
  source(paste(paths$fun, 'TVPRR_VARF_v190304.R', sep='/'))
  
  # Outliers function
  OF = function(pred,y=train[subset,1],tol=2,go.to.pred=pred.lin){
    newx=pred
    
    cond.max = (newx-mean(y))>tol*(max(y)-mean(y))
    cond.min = (newx-mean(y))<tol*(min(y)-mean(y))
    
    newx[cond.max]=go.to.pred[cond.max]
    newx[cond.min]=go.to.pred[cond.min]
    
    return(newx)    
  }
  
  
  # =========================================================================== #
  #                         2. DATA AND POOS PARAMETERS
  # =========================================================================== #
  
  # FRED_MD (stationnarized), dependant variables and date vector
  fred          <- read_csv(paste(paths$dat,'FRED_QD_stationnary.csv',sep='/'))[-1,-1]
  forecast_vars <- read_csv(paste(paths$dat,'newQ_targets.csv',sep='/'))
  forecast_vars <- forecast_vars[-c(1:2),-1] # Drop first two observations
  forecast_vars <- forecast_vars[-c(1:2),] # Drop first two observations
  date          <- seq(as.Date('1959-04-01'), 
                       as.Date('2018-01-01'), 
                       by = "month")
  
  
  
  # Choose variables and horizons and the data set will be created
  vars <- c(6 ,24,  82, 106,92) 
  szv = length(vars)
  
  bigt <- dim(fred)[1]   # Get time dimension
  bign <- dim(fred)[2]   # Get number of variables
  vars = c(1,59,120,147,253) 
  szv = length(vars)
  hor <- c(1,2,4)
  szh <- length(hor) 
  nvar = length(vars) # number of dependent variables to be forecasted
  varsMedium = c(6,24,82,106,92,48,32,2,3,19,100,116,120,65,68,69,64,74,84,95,104) 
  nvar_medium = length(varsMedium)
  factor_haircut <- 60 #how many factors do we consider for LASSO
  
  # Exercise starts
  tau = 158
  bigt=222
  
  # How often to reoptimize hyperparameters
  jump <- 1 # Every two years
  
  # Maximal lags and factors considered
  ly=2
  lf=2
  nf   <- 3
  obstot=360
  silenceplz = 1
  scree_ts=.05
  maxf=3
  alpha=0.15
  sv.param=0
  
  # Sequences
  optim_seq  <- seq(tau,bigt,by=jump)       # Points where we optimize HPs
  compl_seq  <- c((tau):bigt)[-seq(1,length((tau):bigt),by=jump)] # Others
  
  
  var.vars=paste('',c('GDPC1','UNRATE','T10TFFM','CPIAUCSL','HOUST','GS1','PAYEMS','DPIC96','PCECC96','CUMFNS',	'WPSFD49207',
                      'PCECTPI','CES0600000008','M2REAL','TOTRESNS','NONBORRES','M1REAL',	'SP500',	'GS10',	'FEDFUNDS'),sep='')
  vars.pos = c(1,2,4,6,3) 
  
  #get the z.pos
  add.this.all = c() #ar2
  for(j in 1:length(var.vars)){
    add.this=c(grep(var.vars[j], colnames(fred)))
    #double.check=colnames(newtrain)[add.this]==var.vars[j]
    add.this.all=c(add.this.all,add.this)
  }
  varsMedium=add.this.all
  
  
  # =========================================================================== #
  #                                3. POOS
  # =========================================================================== #
  
  # Build dependant variable for inflation models
  M=4
  
  # Initialize matrices for results
  hp_track   <- array(NA, dim=c(bigt,max(hor),5,M,150)) 
  forecast   <- array(NA,dim=c(bigt,max(hor),5,M))
  
  h = all_options$H[it_pos]
  mod = all_options$M[it_pos]
  lambdavec=exp(linspace(-2,12,n=15))
  if(mod==4){lambdavec=lambdavec[4:15]}
  hor_pos = c(1,2,3,3)
  
  sname = paste(paths$out,'/table16to18/TVPfcst_',all_options$V[it_pos],h,mod,'.RData',sep='')
  
  
  for (v in all_options$V[it_pos]){ #all vars
    tic()
    fred2 <- EM_sw(data=fred,n=8,it_max=1000)$data
    
    # START POOS -------------------------------------------------------------
    for(t in tau:bigt){
      
      # INFO_{t-h} (we forecast y_t), excluding related variable for X_{t-h}
      data = as.matrix(fred2)[1:(t-h),]
      y       <- as.matrix(forecast_vars[1:(t-h),(v-1)*5+hor_pos[h]]) # Dependant variable 
      Y       <- as.matrix(forecast_vars[1:(t-h),(v-1)*5+1])     # Regressor          
      factors = data
      
      if(mod==1){
        # Drop missing values
        start   <- sum(is.na(y)) + 1
        end     <- length(1:(t-h))
        y       <- y[start:end,1]      # added second dimension
        Y       <- Y[start:end,1]      # added second dimension
        factors <- factors[start:end,]
        
        # Make regressor matrix and last observations for forecast
        train <- make_reg_matrix(y=y,Y=Y,factors=as.matrix(factors),h=h,ly=ly,lf=lf)[,1:3]
        last = train[dim(train)[1],]
        train = train[1:(dim(train)[1]-h),]
      }
      if(mod==2){
        factors=EM_sw(data=factors[,-vars[v]],n=2,it_max=1000)$factors
        
        # Drop missing values
        start   <- sum(is.na(y)) + 1
        end     <- length(1:(t-h))
        y       <- y[start:end,1]      # added second dimension
        Y       <- Y[start:end,1]      # added second dimension
        factors <- factors[start:end,]
        
        # Make regressor matrix and last observations for forecast
        train <- make_reg_matrix(y=y,Y=Y,factors=as.matrix(factors),h=h,ly=ly,lf=lf)
        last = train[dim(train)[1],]
        train = train[1:(dim(train)[1]-h),]
      }
      if(mod==3){
        factors=factors[,vars[-v]]
        
        # Drop missing values
        start   <- sum(is.na(y)) + 1
        end     <- length(1:(t-h))
        y       <- y[start:end,1]      # added second dimension
        Y       <- Y[start:end,1]      # added second dimension
        factors <- factors[start:end,]
        
        # Make regressor matrix and last observations for forecast
        train <- make_reg_matrix(y=y,Y=Y,factors=as.matrix(factors),h=h,ly=ly,lf=lf)
        last = train[dim(train)[1],]
        train = train[1:(dim(train)[1]-h),]
      }
      if(mod==4){
        
        factors=factors[,varsMedium[-vars.pos[v]]]
        
        # Drop missing values
        start   <- sum(is.na(y)) + 1
        end     <- length(1:(t-h))
        y       <- y[start:end,1]      # added second dimension
        Y       <- Y[start:end,1]      # added second dimension
        factors <- factors[start:end,]
        
        # Make regressor matrix and last observations for forecast
        train <- make_reg_matrix(y=y,Y=Y,factors=as.matrix(factors),h=h,ly=ly,lf=lf)
        last = train[dim(train)[1],]
        train = train[1:(dim(train)[1]-h),]
      }
      
      #Rotation 
      maxlag <- max(lf,ly)
      train=train[(maxlag+1):nrow(train),]
      train= as.data.frame(train)
      train = as.matrix(train[which(complete.cases(train)==TRUE),])
      subset=1:(dim(train)[1]) #-(h-1))
      ts.plot(train[,1])
      print(paste('(h,v,t)=(',h,',',v,',',t,')'))
      
      # =========================================================================== #
      #                 2. ESTIMATE MODELS AND HYPERPARAMETERS
      # =========================================================================== #
      
      # 1. Plain ------------------------------------------------------------------
      CV = cv.glmnet(x=train[subset,-1],y=train[subset,1],  family='gaussian', alpha=0)
      mdl = glmnet(x=train[subset,-1],y=train[subset,1],  family='gaussian', alpha=0,lambda=CV$lambda.min)
      
      m=1
      forecast[t,h,v,m] <- predict(mdl, newx=t(as.matrix(c(last[-1])))) -last[1]
      pred.lin = predict(mdl, newx=t(as.matrix(c(last[-1]))))
      hp_track[t,h,v,m,1]= CV$lambda.min
      
      # 2. 2-step ridge regression (2SRR) -------------------------------------------------
      aa <- TVPRR_cosso(y=train[subset,1],X=train[subset,-1],
                        lambdavec =lambdavec,sweigths=1,type=2,
                        alpha=0.01,silent=silenceplz,kfold=5,lambda2=CV$lambda.min,
                        tol=10^(-6),maxit=10,
                        oosX=last[-1])
      
      m=2
      forecast[t,h,v,m] <- OF(aa$fcast) -last[1]
      hp_track[t,h,v,m,1]= CV$lambda.min 
      hp_track[t,h,v,m,2]= aa$grrats$lambdas[1] 
      starting_betas = aa$grrats$betas_grr #for factor tvp
      starting_betas_2s = aa$grr$betas_grr
      
      # 3. Multistep ridge regression, Sparse TVPs (MSRRs)  ----------------------------------
      aa <- TVPRR(y=train[subset,1],X=train[subset,-1],
                  lambdavec =lambdavec,sweigths=1,type=3,
                  alpha=0.001,silent=silenceplz,kfold=5,lambda2=CV$lambda.min,
                  tol=10^(-5),maxit=15,
                  oosX=last[-1])
      
      m=3
      forecast[t,h,v,m] <- OF(aa$fcast) -last[1]
      hp_track[t,h,v,m,1]= CV$lambda.min 
      hp_track[t,h,v,m,2]= aa$lambda1  
      hp_track[t,h,v,m,3:(2+dim(train)[2]+1-1)]=aa$grra$sigmasq
      
      # 4. Multistep ridge regression, Dense TVPs (MSRRd)  ------------------------------------
      aa <- TVPRR_VARF(Y=train[subset,1],X=train[subset,-1],orthoFac = TRUE,
                       lambdavec =lambdavec,sweigths=1,type=2,fp.model=1,sv.param=sv.param,
                       alpha=alpha,silent=silenceplz,kfold=5,lambda2=CV$lambda.min,max.step.cv = 8,
                       adaptive=1,aparam=-.5,tol=10^(-10),maxit=20,#starting_betas=starting_betas,
                       lambdabooster=1,var.share = scree_ts,override = maxf,id=1,oosX=last[-1])
      
      m=4
      forecast[t,h,v,m] <- OF(aa$fcast) -last[1]
      hp_track[t,h,v,m,]= aa$HPs
      starting_valuesVAR = aa$starter_pack
      
      howlong=toc()
    }
    
    # Save results -----------------------------------------------------------
    save(forecast, hp_track,forecast_vars,
         file=sname)
  }
}


# =============================================================================================================
# FORECASTING
# =============================================================================================================

cl <- makeCluster(ncores)
registerDoParallel(cl)
foreach(it_pos = 1:nrow(all_options)) %dopar% POOS(it_pos)
stopCluster(cl)

