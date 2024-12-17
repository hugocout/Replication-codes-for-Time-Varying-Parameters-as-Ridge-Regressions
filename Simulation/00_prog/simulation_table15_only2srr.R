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

# Array of packages to be checked and installed if not already
myPKGs <- c('bvarsv','pracma','matrixcalc','MASS','glmnet','shrinkTVP','fGarch') 

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
library(fGarch)
library(shrinkTVP)

# Custom Functions -----------------------------------------------

#%# GLM -----------------------------------------------
source(paste(paths$fun, 'evaluator_v181121b.R', sep='/'))

source(paste(paths$too, 'EM_sw.R', sep='/'))
source(paste(paths$too, 'ICp2.R', sep='/'))
source(paste(paths$too, 'Xgenerators_v190127.R', sep='/'))

# Functions for TVPRR
source(paste(paths$fun, 'dualGRRmdA_v190215_v230816.R', sep='/'))
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

# =============================================================================================================
# SIMULATION PARAMETERS
# =============================================================================================================

# All Parameters combinations
combn <- list(O = c(2),
              V = c(1:20),
              H = c(5))
all_options <- expand.grid(combn)
rownames(all_options) <- c()

endofsample<-0*15
kfolds=5
setups <- c(2,1,3,4,5)
crit<-2
models=4
totcol <- 7*8
finalscore <- array(NA,dim=c(3,3,2,4+1+3,totcol))
fracvec = c(.8,.5,0) 
DDvec = c(6,20) 
S=20
distscore <- array(NA,dim=c(3,3,2,4+1,totcol,S))
VIZ=0
scores <- array(NA,dim=c(S,totcol,3,3,4+1+3))

# =============================================================================================================
# SIMULATION
# =============================================================================================================

starttime=Sys.time()

for(it_pos in 1:nrow(all_options)) {
  
  Stype = all_options$H[it_pos]
  if(is.na(Stype)){Stype=1}
  whichsimul = paste('simul',Stype,sep='')
  
  bigT = all_options$O[it_pos]
  if(is.na(bigT)){bigT=2}
  
  tic()
  for(DD in c(1,2,3)[1:2]){
    for(frac in c(3,1,2)[3]){
      for(setup in (setups)){
        
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
            sig_heta=0
          }
        }
        
        for(s in 1:S){
          set.seed(s)
          #simul
          
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
          
          ####### Primiceri, BVAR #########
          evalmat = rep(66,7)
          newstart=2
          
          if(DD==0){
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
          
          ######### ARR FVTP ##########
          lambdavec=exp(linspace(-1,16,n=30))
          XX <-usmacro[1:(T-1),]
          yy <- usmacro[2:T,1] #univariate model
          dimX= M+1 
          
          sweigths = rep(1,dimX)
          sweigths[1]=1
          ######### GRRA ##########
          
          tic()
          
          plot(cv.glmnet(x=XX,y=yy,alpha=0))
          
          aa=tvp.ridge(X=XX,Y=yy,lambda.candidates=lambdavec,oosX=c(),
                       lambda2=cv.glmnet(x=XX,y=yy,alpha=0)$lambda.min,kfold=5,CV.plot=T,CV.2SRR=TRUE,block_size=8,
                       sig.u.param=0.75,sig.eps.param=0.75,ols.prior=1)
          eval = evaluator(truebeta=beta,truesig=truesig,betahat= aa$betas.2srr[1,,],
                           newstart=newstart,endofsample=0,crit=0.0001,visualizar=VIZ)
          plain.lambda1=aa$lambdas 

          source(paste(paths$fun, 'TVPRR_v181111_v230814.R', sep='/'))
          tic()
          
          aa=tvp.ridge(X=XX,Y=yy,lambda.candidates=lambdavec,oosX=c(),
                       lambda2=cv.glmnet(x=scale(XX),y=scale(yy),alpha=0)$lambda.min,kfold=5,CV.plot=T,CV.2SRR=TRUE,block_size=1,
                       sig.u.param=0.75,sig.eps.param=0.75,ols.prior=1)
          eval = evaluator(truebeta=beta,truesig=truesig,betahat= aa$betas.2srr[1,,],
                           newstart=newstart,endofsample=0,crit=0.0001,visualizar=VIZ)
          plain.lambda1=aa$lambdas 

          source(paste(paths$fun, 'TVPRR_v181111_v230814.R', sep='/'))
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
              Kmat_half[begin:end,]=Kmat_half[begin:end,] 
              
              
              if(-nf>-dimX){Kmat_half[begin,]=1000*Kmat_half[begin,]
              }
            }
            Kmat_half[(ncolZ-dimX+1):ncolZ,]=(1/lambdas[2])*Kmat_half[(ncolZ-dimX+1):ncolZ,] 
            
            expmat = (matrix(rexp((ncolZ-dimX)*ncol(Kmat_half)),(ncolZ-dimX),ncol(Kmat_half))) 
            Kmat_half[-c((ncolZ-dimX+1):ncolZ),]=Kmat_half[-c((ncolZ-dimX+1):ncolZ),]
            
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

                  betas_grr[eq,,t] = betas_grr[eq,,t-1]+uhat[positions,eq]
                }
              }
            }
            ############ CI'S ############
            if(CI>0){
              
              Kmh = Kmat_half
              if(param>obs){ #Go Dual
                for(tt in 1:obs){Kmh[,tt]=Kmat_half[,tt]*eweigths[tt]^-1}
              }
              
              if(dimX>nf){
                Kmh=Kmh[1:(dim(Kmat_half)[1]-dimX),]
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
              
              invCov = chol2inv(chol(tcrossprod(Kmh)+diag(dim(Kmh)[1])))
              
              CT= array(1,dim=c(T-1,T-1))
              CT= lower.triangle(CT)
              C=kronecker(diag(nf),CT)
              sandwich=invCov*(sd(y-yhat)^2)*(lambdas[1])^(-1)
              Vb=sqrt(diag(C%*%tcrossprod(sandwich,C))) 
              
              betas_grr_ci <- array(0,dim=c(dim(as.matrix(y))[2],nf,T,2))
              for(c in 1:2){
                for(eq in 1:dim(as.matrix(y))[2]){
                  for(k in 1:nf){ #variables
                    for(t in 1:(T)){
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
          
          
          TVPRR_BBB <- function(B, block.size, x, y, plain.lambda1,uw,oob=FALSE,expo.p=1,l2=0.001,sf=1) {
            OBS=length(y)
            beta.mat <- array(NA, dim = c(B, ncol(x) + 1, OBS))
            W.mat <- array(NA, dim = c(B, OBS))
            y.mat = array(NA, dim = c(B,  OBS))
            
            for (b in 1:B) {
              if (b %in% seq(1, B, 50)) {
                print(b)
              }
              set.seed(b)
              block.size0 = block.size 
              
              w <- (rexp(rate=expo.p, n = round(OBS/ block.size0)))^(-1)
              w <- w[sort(sample(x = c(1:length(w)), size = OBS, replace = TRUE))]
              y.mod <- y 
              x.mod=x 
              
              if(is.null(uw)){
                aa <- TVPRR(y = as.matrix(y.mod), X = as.matrix(x.mod), opt.sige = w^(-1),
                            lambdavec = plain.lambda1 * rexp(1), sweigths = 1,
                            type = 1, silent = 1, kfold = 5, lambda2 = 0.001*rexp(1))
              }
              else{
                lw = 1 
                luw= (rexp(rate=expo.p,ncol(x.mod)+1)) 
                aa <- TVPRR(y = as.matrix(y.mod), X = as.matrix(x.mod), opt.sige = w^(-1), olsprior = 1,
                            lambdavec = plain.lambda1 *lw, 
                            sweigths = uw *luw, 
                            type = 1, silent = 1, kfold = 5, lambda2 = l2*rexp(1)) 
              }
              
              w_trans <- (w) 
              w_trans[w_trans < 0.1] <- 0.1
              W.mat[b,] = w_trans
              
              
              
              beta.mat[b, , ] <- aa$grr$betas_grr[1,, ]
              y.mat[b,]= aa$grr$yhat
              
              if(oob==TRUE){
                #for(kkk in 1:nrow(beta.mat[b, , ])){beta.mat[b,kkk , ] = beta.mat[b,kkk , ] / w_trans}
              }
            }
            
            print(apply(beta.mat[, , 100],2,sd))
            
            percentiles <- array(NA, dim = c(7, ncol(x) + 1, OBS))
            
            for (t in 1:OBS) {
              for (k in 1:(ncol(x) + 1)) {
                #print(W.mat[,t])
                if(oob){
                  # print('c1')
                  percentiles[, k, t] <- DescTools::Quantile(beta.mat[, k, t],weights=1/W.mat[,t], probs = c(0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975)) #,na.rm = TRUE)
                }else{
                  #print('c2')
                  percentiles[, k, t] <- quantile(beta.mat[, k, t], probs = c(0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975),na.rm = TRUE)
                }
              }
            }
            
            print(mean(apply(y.mat,2,sd)))
            
            return(percentiles)
          }
          
          # Usage example
          B <- 150
          bs = 1
          percentiles <- TVPRR_BBB(B, block.size = bs, XX, yy, expo.p= 1,l2=min(cv.glmnet(x=scale(XX),y=scale(yy),alpha=0)$lambda.min,plain.lambda1/length(yy)), 
                                   plain.lambda1,uw=aa$sig.u,oob = FALSE,sf=sd(aa$yhat.2srr[1,]-y)) #,uw)
          

          
          savetoc=toc()
          eval[3]=savetoc
          
          {
            par(mfrow=c(2,3))
            for(vv in 1:nrow(beta)){
              betaa = beta[vv,]
              ts.plot(cbind(betaa,aa$betas.2srr[1,vv,],percentiles[3,vv,],percentiles[5,vv,],percentiles[4,vv,]),col=1:6,lwd=2)
            }
            sin = scale(sig_e[-1])
            if(sd(sig_e)==0){sin=0}

            cove = matrix(NA,3,ncol(XX)+1)
            
            for(vv in 1:(ncol(XX)+1)){
              betaa = beta[vv,]
              cove[1,vv]= mean(I(percentiles[3,vv,]  <betaa & betaa <percentiles[5,vv,]))
              cove[2,vv]= mean(I(percentiles[2,vv,] <betaa & betaa <percentiles[6,vv,]))
              cove[3,vv]= mean(I(percentiles[1,vv,]  <betaa & betaa <percentiles[7,vv,]))
            }
            cov.metrics = apply(cove,1,mean)
            
            print(cov.metrics)
          }
          
          toc()
          
          emp.iqr = apply(percentiles[5,,]-percentiles[3,,],2,mean)
          
          ts.plot(emp.iqr)
          ts.plot(cbind(scale(emp.iqr),sin),col=1:5,lwd=2)
          
          eval = c(eval,cov.metrics,cor(sig_e[-1],emp.iqr))
          evalmat = rbind(evalmat,t(as.vector(eval)))

          ######### COSSO 1 ##########
          eval= rep(66,7)
          evalmat = rbind(evalmat,t(as.vector(eval)))
          
          ######### COSSO F (Oracle) ##########
          eval = rep(66,7)
          
          evalmat = rbind(evalmat,(as.vector(eval)))
          
          
          ############### fancy bayes ###################
          if(DD!=4){

            eval=rep(66,7)
            
            evalmat = rbind(evalmat,(as.vector(eval)))
            
            ############### fancy bayes 2 ###################
            
            data=data.frame(y=yy,XX)
            tic()
            res <- shrinkTVP(y ~ ., data = data,learn_a_xi = FALSE,mod_type = 'ridge', 
                             niter = 20,
                             sv = TRUE)
            savetoc=toc()
            
            #translation
            betas.bayes = matrix(NA,ncol(XX)+1,nrow(XX))
            for(kkk in 1:(ncol(XX)+1)){betas.bayes[kkk,]=apply(res$beta[[kkk]],2,mean)[-(nrow(XX)+1)]}
            
            
            
            eval = evaluator(truebeta=beta,truesig=truesig,betahat= betas.bayes,
                             newstart=newstart,endofsample=0,crit=0.0001,visualizar=VIZ)
            
            eval[3]=savetoc
            
            
            
            
            cove = matrix(NA,3,ncol(XX)+1)
            
            for(vv in 1:(ncol(XX)+1)){
              betaa = beta[vv,]
              kkk=vv
              cove[1,vv]= mean(I(apply(res$beta[[kkk]],2,quantile,probs=0.16)[-(nrow(XX)+1)]  <betaa & betaa < apply(res$beta[[kkk]],2,quantile,probs=0.84)[-(nrow(XX)+1)]))
              cove[2,vv]= mean(I(apply(res$beta[[kkk]],2,quantile,probs=0.05)[-(nrow(XX)+1)]  <betaa & betaa < apply(res$beta[[kkk]],2,quantile,probs=0.95)[-(nrow(XX)+1)]))
              cove[3,vv]= mean(I(apply(res$beta[[kkk]],2,quantile,probs=0.025)[-(nrow(XX)+1)]  <betaa & betaa < apply(res$beta[[kkk]],2,quantile,probs=0.975)[-(nrow(XX)+1)]))
            }
            cov.metrics = apply(cove,1,mean)
            print('bayes')
            print(cov.metrics)
            
            par(mfrow=c(2,3))
            for(vv in 1:nrow(beta)){
              betaa = beta[vv,]
              kkk=vv
            }
            
            emp.iqr = matrix(NA,ncol(XX)+1,nrow(XX))
            for(kkk in 1:(ncol(XX)+1)){
              emp.iqr[kkk,] = apply(res$beta[[kkk]],2,quantile,probs=0.84)[-(nrow(XX)+1)]-apply(res$beta[[kkk]],2,quantile,probs=0.16)[-(nrow(XX)+1)]
            }
            emp.iqr = apply(emp.iqr,2,mean)
            eval = c(eval,cov.metrics,cor(sig_e[-1],emp.iqr))
            
            evalmat = rbind(evalmat,(as.vector(eval)))
            
            ############### fancy bayes 3 ###################
            
            data=data.frame(y=yy,XX)
            tic()
            res <- shrinkTVP(y ~ ., data = data,mod_type = "triple",
                             niter = 20,
                             sv = TRUE)
            savetoc=toc()
            
            #translation
            for(kkk in 1:(ncol(XX)+1)){betas.bayes[kkk,]=apply(res$beta[[kkk]],2,mean)[-(nrow(XX)+1)]}
            
            
            
            eval = evaluator(truebeta=beta,truesig=truesig,betahat= betas.bayes,
                             newstart=newstart,endofsample=0,crit=0.0001,visualizar=VIZ)
            
            
            eval[3]=savetoc
            
            
            
            cove = matrix(NA,3,ncol(XX)+1)
            
            for(vv in 1:(ncol(XX)+1)){
              betaa = beta[vv,]
              kkk=vv
              cove[1,vv]= mean(I(apply(res$beta[[kkk]],2,quantile,probs=0.16)[-(nrow(XX)+1)]  <betaa & betaa < apply(res$beta[[kkk]],2,quantile,probs=0.84)[-(nrow(XX)+1)]))
              cove[2,vv]= mean(I(apply(res$beta[[kkk]],2,quantile,probs=0.05)[-(nrow(XX)+1)]  <betaa & betaa < apply(res$beta[[kkk]],2,quantile,probs=0.95)[-(nrow(XX)+1)]))
              cove[3,vv]= mean(I(apply(res$beta[[kkk]],2,quantile,probs=0.025)[-(nrow(XX)+1)]  <betaa & betaa < apply(res$beta[[kkk]],2,quantile,probs=0.975)[-(nrow(XX)+1)]))
            }
            cov.metrics = apply(cove,1,mean)
            print('bayes')
            print(cov.metrics)
            
            betaa = beta[2,]
            kkk=2
     
            emp.iqr = matrix(NA,ncol(XX)+1,nrow(XX))
            for(kkk in 1:(ncol(XX)+1)){
              emp.iqr[kkk,] = apply(res$beta[[kkk]],2,quantile,probs=0.84)[-(nrow(XX)+1)]-apply(res$beta[[kkk]],2,quantile,probs=0.16)[-(nrow(XX)+1)]
            }
            emp.iqr = apply(emp.iqr,2,mean)
            eval = c(eval,cov.metrics,cor(sig_e[-1],emp.iqr))
            evalmat = rbind(evalmat,(as.vector(eval)))
            
          }else{
            evalmat = rbind(evalmat,rep(66,7))
            evalmat = rbind(evalmat,rep(66,7))
            evalmat = rbind(evalmat,rep(66,7))
          }
          
          ########### BB 2SRR ############
          eval=rep(66,7)
          evalmat = rbind(evalmat,t(as.vector(eval)))
          
          ######### Comparaison ##########
          scores[s,1:totcol,DD,frac,setup]=t((vec(evalmat)))
          print('ahhhh')
          print(mean(scores[1:s,26,DD,frac,setup]))
        }}
      
      print(evalmat)
      filename= paste(paths$out,'/table15/output_',whichsimul,'_230817bands_2ssr.Rdata',sep='')
      save(finalscore,scores,distscore,file=filename)
      
    }
  }
  toc()
}

print(Sys.time()-starttime)

