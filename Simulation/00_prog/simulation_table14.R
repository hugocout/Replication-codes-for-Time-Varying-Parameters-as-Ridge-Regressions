# =============================================================================================================
# INITIALIZATION
# =============================================================================================================

rm(list=ls()) # To clear the current environment
yymmdd = paste(substr(Sys.Date(),3,4),substr(Sys.Date(),6,7),substr(Sys.Date(),9,10),sep='')

wd = 'INSERT YOUR PATH HERE/Simulation/'
setwd(wd)
paths <- list(pro = "00_prog",
              too = "10_tools",
              fun = "10_tools/functions",
              sim = "10_tools/simul_types",
              out = "20_output",
              rst = "30_results")

# Functions & Libraries
source(paste(paths$fun, 'evaluator_v181121b.R', sep='/'))
source(paste(paths$too, 'EM_sw.R', sep='/'))
source(paste(paths$too, 'ICp2.R', sep='/'))
source(paste(paths$too, 'Xgenerators_v190127.R', sep='/'))

# Functions for TVPRR
source(paste(paths$fun, 'dualGRRmdA_v190215.R', sep='/'))
source(paste(paths$fun, 'CVGSBHK_v181127.R', sep='/'))
source(paste(paths$fun, 'CVKFMV_v190214.R', sep='/'))
source(paste(paths$fun, 'zfun_v190304.R', sep='/'))
source(paste(paths$fun, 'factor.R', sep='/'))
source(paste(paths$fun, 'TVPRRcosso_v181120.R', sep='/'))
source(paste(paths$fun, 'TVPRRcossoF_v190125.R', sep='/'))
source(paste(paths$fun, 'TVPRR_v181111.R', sep='/'))
source(paste(paths$fun, 'fastZrot_v181125.R', sep='/'))
source(paste(paths$fun, 'TVPRR_VARF_v190304.R', sep='/'))
source(paste(paths$fun, 'MV2SRR_v221103_forCC.R', sep='/'))

#Install
# Array of packages to be checked and installed if not already
myPKGs <- c('bvarsv','pracma','matrixcalc','MASS','glmnet','timeSeries','fGarch','shrinkTVP')

# check to see if each package is installed, and if not add it to a list to install
InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]

# install any needed packages
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

library(bvarsv)
library(pracma)
library(matrixcalc)
library(MASS)
library(glmnet)
library(timeSeries)
library(fGarch)
library(shrinkTVP)


# =============================================================================================================
# SIMULATION : Estimating a Single Equation
# =============================================================================================================

# Parameters ----------------------------------------------------------------
Stype = 5
totcol <- 3*8
finalscore <- array(NA,dim=c(3,3,2,4+1+3,totcol))
S=5
scores <- array(NA,dim=c(S,totcol,4,3,4+1+3))

# Simulation start ----------------------------------------------------------
for(bigT in c(2,1,3)){
  
  
  whichsimul = paste('simul',Stype,sep='')
  endofsample<-0*15
  kfolds=5
  setups <- c(2,1,3,4,5)[1] 
  crit<-2
  models=4
  fracvec = c(.8,.5,0) 
  DDvec = c(6,20,50,100) 
  distscore <- array(NA,dim=c(3,3,2,4+1,totcol,S))
  VIZ=0
  
  starttime=Sys.time()
  
  for(DD in c(1,2,3,4)){
    for(frac in c(3,1,2)[2]){
      for(setup in setups){
        
        obs <- c(151,301,601)[bigT]
        T=obs
        L<-1
        Mbase<-DDvec[DD]
        compressor = .5
        #setups
        {
          if( setup ==1){
            M <-Mbase
            sig_e = 2*rep(1,obs)
            sig_heta=0
          }
          if( setup ==2){
            M <-Mbase
            sig_e = 4*rep(1,obs)
            sig_heta=0
          }
          if( setup ==3){
            M <-Mbase
            sig_e = 6*1.5*rep(1,obs)
            sig_heta=0
          }
          if( setup ==4){
            M <-Mbase
            vv <- as.matrix(rnorm(obs))/10
            sig_e = 0.1*rep(1,obs)
            sig_e[1]=0.5+ (vv[1])
            for(t in 2:obs){sig_e[t]=0.99*sig_e[t-1]+vv[t]}
            sig_e = 0.1*exp(sig_e)
            sig_e = 3*sig_e/mean(sig_e) #so it has mean 1
            sig_e = (sig_e-min(sig_e))/(max(sig_e)-min(sig_e))
            sig_e = 2+(4-2)*sig_e
            #ts.plot(sig_e)
            sig_heta=0
          }
          if( setup ==5){
            M <-Mbase
            vv <- as.matrix(rnorm(obs))/10
            sig_e = 0.1*rep(1,obs)
            sig_e[1]=0.5+ (vv[1])
            for(t in 2:obs){sig_e[t]=0.99*sig_e[t-1]+vv[t]}
            sig_e = 0.1*exp(sig_e)
            sig_e = 3*sig_e/mean(sig_e) #so it has mean 1
            sig_e = (sig_e-min(sig_e))/(max(sig_e)-min(sig_e))
            sig_e = 2+(1.5*6-2)*sig_e
            #ts.plot(sig_e)
            sig_heta=0
          }
        }
        
        for(s in 1:S){
          set.seed(s)
          
          rule= 1
          sori=s
          iii=0
          while(rule>0.95){
            source(paste(paths$sim,'/',whichsimul,'.R',sep=''))
            rule = cor(yhat,y)^2
            iii=1
            if(rule>0.95){s=s+100*iii}
            print(rule)
          }
          s=sori
          
          #TRACKER
          print(paste('We are at (S,D,frac,setup,s)=(',Stype,DDvec[DD],',',fracvec[frac],',',setup,',',s,')'))
          print(cor(yhat,y)^2)
          #}}}}
          
          # Empty eval -------------------------------------------------------------
          evalmat = rep(66,3)
          newstart=2
          
          
          # ARR FVTP ----------------------------------------------------------------
          lambdavec=exp(linspace(2,16,n=15))
          XX <-usmacro[1:(T-1),]
          yy <- usmacro[2:T,1] #univariate model
          dimX= M+1 # + constant
          
          sweigths = rep(1,dimX)
          sweigths[1]=1
          
          # 2SRR --------------------------------------------------------------------
          
          tic()
          aa=tvp.ridge(X=XX,Y=yy,lambda.candidates=lambdavec,oosX=c(),
                       lambda2=.0001,kfold=5,CV.plot=F,CV.2SRR=TRUE,block_size=1,
                       sig.u.param=.75,sig.eps.param=.75,ols.prior=0)
          savetoc=toc()
          
          eval = evaluator(truebeta=beta,truesig=truesig,betahat= aa$betas.2srr[1,,],
                           newstart=newstart,endofsample=0,crit=0.0001,visualizar=VIZ)
          
          eval[3]=savetoc
          evalmat = rbind(evalmat,t(as.vector(eval)))
          
          # Empty eval -------------------------------------------------------------
          eval = rep(66,3)
          evalmat = rbind(evalmat,(as.vector(eval)))
          
          # Empty eval -------------------------------------------------------------
          eval = rep(66,3)
          evalmat = rbind(evalmat,(as.vector(eval)))
          
          # Empty eval -------------------------------------------------------------
          evalmat = rbind(evalmat,rep(66,3))
          evalmat = rbind(evalmat,rep(66,3))
          evalmat = rbind(evalmat,rep(66,3))
          
          
          # Empty eval -------------------------------------------------------------
          evalmat = rbind(evalmat,rep(66,3))
          
          # Comparaison -------------------------------------------------------------
          scores[s,1:totcol,DD,bigT,setup]=t((vec(evalmat)))
          print(evalmat)
        }
      }
    }
  }
}
print(Sys.time()-starttime)


# Results -------------------------------------------------------------------
results_toppanel = apply(scores[,c(18),,,2],c(2,3),mean,na.rm=TRUE)
colnames(results_toppanel) = c("150","300","600")
rownames(results_toppanel) = c("Single equation, 6","20","50","100")

# =============================================================================================================
# SIMULATION : Estimating K Equations
# =============================================================================================================

# Parameters ----------------------------------------------------------------
Stype = 5
totcol <- 3*8
finalscore <- array(NA,dim=c(3,3,2,4+1+3,totcol))
S=5
scores <- array(NA,dim=c(S,totcol,4,3,4+1+3))

# Simulation start ----------------------------------------------------------
for(bigT in c(2,1,3)){
  
  
  whichsimul = paste('simul',Stype,sep='')
  endofsample<-0*15
  kfolds=5
  setups <- c(2,1,3,4,5)[1] 
  crit<-2
  models=4
  fracvec = c(.8,.5,0) 
  DDvec = c(6,20,50,100) 
  distscore <- array(NA,dim=c(3,3,2,4+1,totcol,S))
  VIZ=0
  
  starttime=Sys.time()
  
  for(DD in c(1,2,3,4)){
    for(frac in c(3,1,2)[2]){
      for(setup in setups){
        #print(setup)
        obs <- c(151,301,601)[bigT]
        T=obs
        L<-1
        Mbase<-DDvec[DD]
        compressor = .5
        #setups
        {
          if( setup ==1){
            M <-Mbase
            sig_e = 2*rep(1,obs)
            sig_heta=0
          }
          if( setup ==2){
            M <-Mbase
            sig_e = 4*rep(1,obs)
            sig_heta=0
          }
          if( setup ==3){
            M <-Mbase
            sig_e = 6*1.5*rep(1,obs)
            sig_heta=0
          }
          if( setup ==4){
            M <-Mbase
            vv <- as.matrix(rnorm(obs))/10
            sig_e = 0.1*rep(1,obs)
            sig_e[1]=0.5+ (vv[1])
            for(t in 2:obs){sig_e[t]=0.99*sig_e[t-1]+vv[t]}
            sig_e = 0.1*exp(sig_e)
            sig_e = 3*sig_e/mean(sig_e) #so it has mean 1
            sig_e = (sig_e-min(sig_e))/(max(sig_e)-min(sig_e))
            sig_e = 2+(4-2)*sig_e
            #ts.plot(sig_e)
            sig_heta=0
          }
          if( setup ==5){
            M <-Mbase
            vv <- as.matrix(rnorm(obs))/10
            sig_e = 0.1*rep(1,obs)
            sig_e[1]=0.5+ (vv[1])
            for(t in 2:obs){sig_e[t]=0.99*sig_e[t-1]+vv[t]}
            sig_e = 0.1*exp(sig_e)
            sig_e = 3*sig_e/mean(sig_e) #so it has mean 1
            sig_e = (sig_e-min(sig_e))/(max(sig_e)-min(sig_e))
            sig_e = 2+(1.5*6-2)*sig_e
            #ts.plot(sig_e)
            sig_heta=0
          }
        }
        
        for(s in 1:S){
          if(is.na(s)){s=1}
          set.seed(s)
          
          # simul
          rule= 1
          sori=s
          iii=0
          while(rule>0.95){
            source(paste(paths$sims,'/',whichsimul,'.R',sep=''))
            rule = cor(yhat,y)^2
            iii=1
            if(rule>0.95){s=s+100*iii}
            print(rule)
          }
          s=sori
          
          #TRACKER
          print(paste('We are at (S,D,frac,setup,s)=(',Stype,DDvec[DD],',',fracvec[frac],',',setup,',',s,')'))
          print(cor(yhat,y)^2)
          #}}}}
          
          # Primiceri, BVAR ------------------------------------------------
          evalmat = rep(66,3)
          newstart=2
          if(DD==666){
            tau=50
            newstart=tau+1
            tic()
            bvar_est <- bvar.sv.tvp(usmacro, p = 1, tau = tau, nf = 10, pdrift = FALSE, nrep = 10000,
                                    nburn = 5000, thinfac = 10, itprint = 500000, save.parameters = TRUE,
                                    k_B = 4, k_A = 4, k_sig = 1, k_Q = 0.01, k_S = 0.1, k_W = 0.01,
                                    pQ = NULL, pW = NULL, pS = NULL)
            savetoc=toc()
            stop
            eval = evaluator(truebeta=beta[,newstart:dim(beta)[2]],truesig=truesig,betahat= bvar_est$Beta.postmean[1,,],
                             newstart=2,endofsample=0,crit=0.0001,visualizar=VIZ)
            eval[3]=savetoc
            evalmat= t(as.vector(eval))
          }
          
          # ARR FVTP ------------------------------------------------------
          lambdavec=exp(linspace(2,16,n=15))
          XX <-usmacro[1:(T-1),]
          yy <- usmacro[2:T,] #univariate model
          dimX= M+1 # + constant
          
          sweigths = rep(1,dimX)
          sweigths[1]=1
          
          # 2SRR ----------------------------------------------------------
          
          tic()
          aa=tvp.ridge(X=XX,Y=yy,lambda.candidates=lambdavec,oosX=c(),
                       lambda2=.0001,kfold=5,CV.plot=F,CV.2SRR=FALSE,block_size=1,
                       sig.u.param=.75,sig.eps.param=.75,ols.prior=0)
          savetoc=toc()
          
          
          
          eval = evaluator(truebeta=beta,truesig=truesig,betahat= aa$betas.2srr[1,,],
                           newstart=newstart,endofsample=0,crit=0.0001,visualizar=VIZ)
          
          eval[3]=savetoc
          evalmat = rbind(evalmat,t(as.vector(eval)))
          
          ##
          if(DD==666){
            tic()
            aa <- TVPRR(y=yy,X=XX,
                        lambdavec =lambdavec,sweigths=sweigths,type=3,
                        alpha=0.5,silent=abs(VIZ-1),kfold=5,lambda2=0.0001)
            savetoc=toc()
            
            
            eval = evaluator(truebeta=beta,truesig=truesig,betahat= aa$grra$betas_grr[1,,],
                             newstart=newstart,endofsample=0,crit=0.0001,visualizar=VIZ)
            eval[3]=savetoc
            
            evalmat = rbind(evalmat,t(as.vector(eval)))
          }
          
          eval = rep(66,3)
          
          evalmat = rbind(evalmat,(as.vector(eval)))
          
          # Empty eval -------------------------------------------------------------
          eval = rep(66,3)
          
          evalmat = rbind(evalmat,(as.vector(eval)))
          
          
          # Empty eval -------------------------------------------------------------
          evalmat = rbind(evalmat,rep(66,3))
          evalmat = rbind(evalmat,rep(66,3))
          evalmat = rbind(evalmat,rep(66,3))
          
          
          # Empty eval -------------------------------------------------------------
          evalmat = rbind(evalmat,rep(66,3))
          
          # Comparaison -------------------------------------------------------------
          scores[s,1:totcol,DD,bigT,setup]=t((vec(evalmat)))
          print(evalmat)
        }
        
      }
    }
  }
}

print(Sys.time()-starttime)

# Results -------------------------------------------------------------------
results_bottompanel = apply(scores[,c(18),,,setup],c(2,3),mean,na.rm=TRUE)
colnames(results_bottompanel) = c("150","300","600")
rownames(results_bottompanel) = c("K equations, 6","20","50","100")
print(results_bottompanel)

# =============================================================================================================
# RESULTS
# =============================================================================================================

table = rbind(results_toppanel,results_bottompanel)
print(table)

# Save to csv
write.csv(table, paste(paths$rst, '/table14.csv', sep='/'), row.names = TRUE)
