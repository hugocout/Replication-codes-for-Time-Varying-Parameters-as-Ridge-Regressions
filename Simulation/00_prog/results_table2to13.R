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

# Array of packages to be checked and installed if not already
myPKGs <- c('pracma','Hmisc')

# check to see if each package is installed, and if not add it to a list to install
InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]

# install any needed packages
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

# Load packages
library(pracma)
library(Hmisc)

# =============================================================================================================
# LOAD RESULTS : Tables 2 to 5
# =============================================================================================================

table_names <- c("table2","table3","table4","table5")
postable <- 1
for(whatS in c(0,3,4,5)){
  
  # Create storage
  load(paste0(paths$out,'/table2to13/output_simul',whatS,'_210308_300obs_s2.Rdata'))
  DIM=c(dim(scores))
  DIM[1]=100
  score2 = array(NA,dim=DIM)
  
  # Load all results
  for(s in c(1:50)){
    load(paste0(paths$out,'/table2to13/output_simul',whatS,'_210308_300obs_s',s,'.Rdata'))
    chosen=scores[s,,,,]
    print(length(chosen[(chosen>10 & chosen!=66)]))
    score2[s,,,,]=chosen
  } 
  scores=score2
  rm(score)
  totcol=15
  finalscore <- array(NA,dim=c(3,3,2,5,totcol))
  
  # Compute mean and median
  for(DD in c(1,2,3)){
    for(frac in c(1:3)){
      for(setup in 1:5){
        #print(scores[,,DD,frac,setup])
        for(j in 1:totcol){
          finalscore[DD,frac,1,setup,j] = mean(scores[,j,DD,frac,setup],na.rm=TRUE)
          finalscore[DD,frac,2,setup,j] = median(scores[,j,DD,frac,setup],na.rm = TRUE)
        }
      }}}
  
  DD=3
  frac=2
  setup=3
  
  # Make Latex table ------------------------------------------------------------------
  
  fracs = c('$\\mathbf{\\sfrac{K^*}{K}=0.2}$','$\\mathbf{\\sfrac{K^*}{K}=0.5}$','$\\mathbf{\\sfrac{K^*}{K}=1}$')
  sigmas = c('$\\sigma_{\\epsilon} =\\text{Low}$','$\\sigma_{\\epsilon} =\\text{Medium}$','$\\sigma_{\\epsilon} =\\text{High}$',
             '$\\sigma_{\\epsilon,t} = \\text{SV Low-Med}$','$\\sigma_{\\epsilon,t}  = \\text{SV Low-High}$')
  TVnames = c('BVAR','2SRR','MSRRs','MSRRd')
  Kvec = c('$K=6$','$K=20$','$K=100$')
  
  
  # Simulation 1, no corr,alpha=.5
  newtmp =finalscore[,,1,,c(1,3:5)]
  newtmp[newtmp>10] = 66
  newtmp = format(round(newtmp,digits=3),nsmall=3)
  newtmp[is.na(newtmp)]='-'
  newtmp[newtmp=='66.000'] = '-'
  
  ft <- ftable(newtmp,row.vars = c(2,3),col.vars = c(1,4)) 
  print(ft)
  
  newtmp[newtmp==66] = '-'
  dimnames(newtmp)[[1]]=Kvec
  dimnames(newtmp)[[2]]=fracs
  dimnames(newtmp)[[3]]=sigmas
  dimnames(newtmp)[[4]]=TVnames
  newtmp[is.na(newtmp)]='-'
  ft <- ftable(newtmp,row.vars = c(2,3),col.vars = c(1,4)) 
  ft
  
  filename <- paste0(paths$rst,"/",table_names[postable],'.tex')
  latex(ft, caption=paste('Results for Simulation ',whatS,sep=''),
        colheads=rep(attr(ft,"col.vars")[[2]],4), 
        rowname=attr(ft,"row.vars")[[2]], label='s1_table',
        rgroup=attr(ft,"row.vars")[[1]],
        cgroup=attr(ft,"col.vars")[[1]],
        title='',file=filename)
  
  postable = postable+1
  
}

# =============================================================================================================
# LOAD RESULTS : Tables 6 to 9
# =============================================================================================================

table_names <- c("table6","table7","table8","table9")
postable <- 1
for(whatS in c(0,3,4,5)){
  
  # Create storage
  load(paste0(paths$out,'/table2to13/output_simul',whatS,'_210308_150obs_s2.Rdata'))
  DIM=c(dim(scores))
  DIM[1]=100
  score2 = array(NA,dim=DIM)
  
  # Load all results
  for(s in c(1:50)){
    load(paste0(paths$out,'/table2to13/output_simul',whatS,'_210308_150obs_s',s,'.Rdata'))
    chosen=scores[s,,,,]
    print(length(chosen[(chosen>10 & chosen!=66)]))
    score2[s,,,,]=chosen
  } 
  scores=score2
  rm(score)
  totcol=15
  finalscore <- array(NA,dim=c(3,3,2,5,totcol))
  
  # Compute mean and median
  for(DD in c(1,2,3)){
    for(frac in c(1:3)){
      for(setup in 1:5){
        #print(scores[,,DD,frac,setup])
        for(j in 1:totcol){
          finalscore[DD,frac,1,setup,j] = mean(scores[,j,DD,frac,setup],na.rm=TRUE)
          finalscore[DD,frac,2,setup,j] = median(scores[,j,DD,frac,setup],na.rm = TRUE)
        }
      }}}
  
  DD=3
  frac=2
  setup=3
  
  # Make Latex table ------------------------------------------------------------------
  
  fracs = c('$\\mathbf{\\sfrac{K^*}{K}=0.2}$','$\\mathbf{\\sfrac{K^*}{K}=0.5}$','$\\mathbf{\\sfrac{K^*}{K}=1}$')
  sigmas = c('$\\sigma_{\\epsilon} =\\text{Low}$','$\\sigma_{\\epsilon} =\\text{Medium}$','$\\sigma_{\\epsilon} =\\text{High}$',
             '$\\sigma_{\\epsilon,t} = \\text{SV Low-Med}$','$\\sigma_{\\epsilon,t}  = \\text{SV Low-High}$')
  TVnames = c('BVAR','2SRR','MSRRs','MSRRd')
  Kvec = c('$K=6$','$K=20$','$K=100$')
  
  
  # Simulation 1, no corr,alpha=.5
  newtmp =finalscore[,,1,,c(1,3:5)]
  newtmp[newtmp>10] = 66
  newtmp = format(round(newtmp,digits=3),nsmall=3)
  newtmp[is.na(newtmp)]='-'
  newtmp[newtmp=='66.000'] = '-'
  
  ft <- ftable(newtmp,row.vars = c(2,3),col.vars = c(1,4)) 
  print(ft)
  
  newtmp[newtmp==66] = '-'
  dimnames(newtmp)[[1]]=Kvec
  dimnames(newtmp)[[2]]=fracs
  dimnames(newtmp)[[3]]=sigmas
  dimnames(newtmp)[[4]]=TVnames
  newtmp[is.na(newtmp)]='-'
  ft <- ftable(newtmp,row.vars = c(2,3),col.vars = c(1,4)) 
  ft
  
  filename <- paste0(paths$rst,"/",table_names[postable],'.tex')
  latex(ft, caption=paste('Results for Simulation ',whatS,sep=''),
        colheads=rep(attr(ft,"col.vars")[[2]],4), 
        rowname=attr(ft,"row.vars")[[2]], label='s1_table',
        rgroup=attr(ft,"row.vars")[[1]],
        cgroup=attr(ft,"col.vars")[[1]],
        title='',file=filename)
  
  postable = postable+1
  
}

# =============================================================================================================
# LOAD RESULTS : Tables 10 to 13
# =============================================================================================================

table_names <- c("table10","table11","table12","table13")
postable <- 1
for(whatS in c(0,3,4,5)){
  
  # Create storage
  load(paste0(paths$out,'/table2to13/output_simul',whatS,'_210308_150obs_s2.Rdata'))
  DIM=c(dim(scores))
  DIM[1]=100
  score2 = array(NA,dim=DIM)
  
  # Load all results
  for(s in c(1:50)){
    load(paste0(paths$out,'/table2to13/output_simul',whatS,'_210308_150obs_s',s,'.Rdata'))
    chosen=scores[s,,,,]
    print(length(chosen[(chosen>10 & chosen!=66)]))
    score2[s,,,,]=chosen
  } 
  scores=score2
  rm(score)
  totcol=15
  finalscore <- array(NA,dim=c(3,3,2,5,totcol))
  
  # Compute mean and median
  for(DD in c(1,2,3)){
    for(frac in c(1:3)){
      for(setup in 1:5){
        #print(scores[,,DD,frac,setup])
        for(j in 1:totcol){
          finalscore[DD,frac,1,setup,j] = mean(scores[,j,DD,frac,setup],na.rm=TRUE)
          finalscore[DD,frac,2,setup,j] = median(scores[,j,DD,frac,setup],na.rm = TRUE)
        }
      }}}
  
  DD=3
  frac=2
  setup=3
  
  # Make Latex table ------------------------------------------------------------------
  
  fracs = c('$\\mathbf{\\sfrac{K^*}{K}=0.2}$','$\\mathbf{\\sfrac{K^*}{K}=0.5}$','$\\mathbf{\\sfrac{K^*}{K}=1}$')
  sigmas = c('$\\sigma_{\\epsilon} =\\text{Low}$','$\\sigma_{\\epsilon} =\\text{Medium}$','$\\sigma_{\\epsilon} =\\text{High}$',
             '$\\sigma_{\\epsilon,t} = \\text{SV Low-Med}$','$\\sigma_{\\epsilon,t}  = \\text{SV Low-High}$')
  TVnames = c('BVAR','2SRR','MSRRs','MSRRd')
  Kvec = c('$K=6$','$K=20$','$K=100$')
  
  
  # Simulation 1, no corr,alpha=.5
  newtmp =finalscore[,,1,,c(1,3:5)]
  newtmp[newtmp>10] = 66
  newtmp = format(round(newtmp,digits=3),nsmall=3)
  newtmp[is.na(newtmp)]='-'
  newtmp[newtmp=='66.000'] = '-'
  
  ft <- ftable(newtmp,row.vars = c(2,3),col.vars = c(1,4)) 
  print(ft)
  
  newtmp[newtmp==66] = '-'
  dimnames(newtmp)[[1]]=Kvec
  dimnames(newtmp)[[2]]=fracs
  dimnames(newtmp)[[3]]=sigmas
  dimnames(newtmp)[[4]]=TVnames
  newtmp[is.na(newtmp)]='-'
  ft <- ftable(newtmp,row.vars = c(2,3),col.vars = c(1,4)) 
  ft
  
  filename <- paste0(paths$rst,"/",table_names[postable],'.tex')
  latex(ft, caption=paste('Results for Simulation ',whatS,sep=''),
        colheads=rep(attr(ft,"col.vars")[[2]],4), 
        rowname=attr(ft,"row.vars")[[2]], label='s1_table',
        rgroup=attr(ft,"row.vars")[[1]],
        cgroup=attr(ft,"col.vars")[[1]],
        title='',file=filename)
  
  postable = postable+1
  
}
