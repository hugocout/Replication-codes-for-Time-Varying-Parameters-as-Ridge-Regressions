
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
myPKGs <- c('pracma','Hmisc','forecast')

# check to see if each package is installed, and if not add it to a list to install
InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]

# install any needed packages
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

# Libraries
library(pracma)     # For repmat 
library(Hmisc)
library(forecast)   # Diebold-Mariano test

# Remove scientific notation
options(scipen = 999)

# =============================================================================================================
# LOAD RESULTS
# =============================================================================================================

# Load and create storage
load(paste(paths$out,'/table16to18/TVPfcst_111_block.Rdata',sep=''))
dim_forecast <- dim(forecast)
dim_forecast[2] <- 4
old_err = array(NA,dim=c(dim_forecast,4))

for(mod in 1:4){
  for(v in 1:5){
    for(h in c(1,2,4)){
      
      # Load
      load(paste(paths$out,'/table16to18/TVPfcst_',v,h,mod,'_block.Rdata',sep=''))
      h_pos <- h
      if(h == 4) {
        h_pos <- dim(forecast)[2]
      }
      
      # Store table results
      old_err[,h,v,,mod]=forecast[,h_pos,v,]
      
    }
  }
}

old_err=old_err[158:222,-3,,,] # Remove empty columns

# =============================================================================================================
# DIEBOLD-MARIANO TESTS (vs AR base model)
# =============================================================================================================

# Parameters --------------------------------------------------------------------------------
M=18
V=5
H=3
compare.to=1
hor=c(1,2,4)

# D-M Test for each models -------------------------------------------------------------------
DMdo=1
if(DMdo ==1){
  DM_mat <- array(0,dim=c(H,V,M,4))
  minM=1
  
  for(v in 1:V){
    for(h in 1:H){
      for(m in 1:4){
        for(mod  in 1:4){
          if(!(m==1 & mod==1)){
            DM <- dm.test( old_err[,h,v,m,mod],old_err[,h,v,1,1], alternative = "two.sided", h = hor[h], power = 2)
            DM_mat[h,v,m,mod] <- DM$p.value
          }
          else{DM_mat[h,v,m,mod] <- 1}
        }
      }
    }
  } 
}

# =============================================================================================================
# TABLE PROTOTYPE
# =============================================================================================================

# Find the overall best model for a (v,h) combinaison -------------------------------------------
mspe_mattemp = sqrt(apply(old_err[,,,,]^2,c(2,3,4,5),mean)) 
mspe_mat = mspe_mattemp
mins = array(0,dim=c(H,V,2))

for(h in 1:H){
  for(v in 1:V){
    mins[h,v,] = which(mspe_mattemp[h,v,,] == min(mspe_mattemp[h,v,,]), arr.ind = TRUE)[1,]
  }
}

for(m in 1:4){
  for(mod in 1:4){
    mspe_mat[,,m,mod] = mspe_mattemp[,,m,mod]/mspe_mattemp[,,1,1]
  }
} 

# Store results
mspe_MAT = format(round(mspe_mat[,,,],digits=2),nsmall=2) 
fcstres=mspe_MAT
MINS=mins

# Attemp at cond tables
for(h in 1:H){for(v in 1:V){for(m in 1:4){for(mod in 1:4){
  fcstres[h,v,m,mod]=mspe_MAT[h,v,m,mod]
  if(-DM_mat[h,v,m,mod]>-0.1){fcstres[h,v,m,mod]=paste(mspe_MAT[h,v,m,mod],'*',sep='')}
  if(-DM_mat[h,v,m,mod]>-0.05){fcstres[h,v,m,mod]=paste(mspe_MAT[h,v,m,mod],'**',sep='')}
  if(-DM_mat[h,v,m,mod]>-0.01){fcstres[h,v,m,mod]=paste(mspe_MAT[h,v,m,mod],'***',sep='')}
}}}}

# Make results table -----------------------------------------------------------------------
fff = fcstres
dimnames(fcstres)[[1]]=c('$h=1$','$h=2$','$h=4$')
dimnames(fcstres)[[2]]=c('GDP','UR','INF','IR','SPREAD')
dimnames(fcstres)[[3]]=c('Plain','2SRR','MSRRs','MSRRd')
dimnames(fcstres)[[4]]=c('AR','ARDI','VAR5','VAR20')

ft_fcst <- ftable(fcstres,row.vars = c(2,1),col.vars = c(4,3)) 
#print(ft_fcst)
ft_fcst1=ft_fcst

# Save to Latex -----------------------------------------------------------------------
extension='.tex'
filename= paste(paths$rst, '/RMSE_table',18,extension,sep='')
ft = ft_fcst1
dviname <- latex(ft, landscape = TRUE, caption='',
                 colheads=rep(attr(ft,"col.vars")[[2]],4),
                 rowname=attr(ft,"row.vars")[[2]], 
                 rgroup=attr(ft,"row.vars")[[1]], 
                 cgroup=attr(ft,"col.vars")[[1]],
                 title='',file=filename) 

# =============================================================================================================
# PRINT RESULTS
# =============================================================================================================

cat("Print Table 18 results:")
print(ft_fcst1)

