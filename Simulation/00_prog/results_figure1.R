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
myPKGs <- c('pracma','Hmisc','RColorBrewer')

# check to see if each package is installed, and if not add it to a list to install
InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]

# install any needed packages
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

# Load packages
library(pracma)
library(Hmisc)
library(RColorBrewer)

# =============================================================================================================
# LOAD RESULTS : Tables 2 to 5
# =============================================================================================================

finalscore2 <- array(NA,dim=c(100,3,3,5,4))

iii=0
for(whatS in c(0,3,4,5)){
  iii=iii+1
  
  load(paste(paths$out,'/table2to13/output_simul',whatS,'_210308_300obs_s27.Rdata',sep=''))
  score = array(NA,dim=c(dim(scores)))
  

  for(s in 1:50){
    
    # Trycatch to avoid errors
    tryCatch({
      load(paste(paths$out,'/table2to13/output_simul',whatS,'_210308_300obs_s',s,'.Rdata',sep=''))
      chosen=scores[s,,,,]
      print(length(chosen[(chosen>10 & chosen!=66)]))
      chosen[(chosen>10 & chosen!=66)]=NA
      score[s,,,,]=chosen
    }, error = function(e) {
      print(paste('Error in simulation',s))
    })
  } 
  scores=score
  rm(score)
  totcol=15
  
  
  for(DD in c(1,2,3)){
    for(frac in c(1:3)){
      for(setup in 1:5){
        #print(scores[,,DD,frac,setup])
        for(j in 3){
          finalscore2[,DD,frac,setup,iii] = scores[,j,DD,frac,setup]/scores[,1,DD,frac,setup]
        }
      }}}
}

finalscore=finalscore2

# =============================================================================================================
# MAKE FIGURE 1 : Boxplots
# =============================================================================================================

# Plot parameters
cols=brewer.pal(n = 7, name = 'Dark2')
par(mfrow=c(1,1), mai = 2*c(0.2, 0.3, 0.1, 0.1))
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))


# Plot (a) By noise process -------------------------------------------------------------------------

# Data
two.vectors = cbind(as.vector(finalscore[,1,,1,]),as.vector(finalscore[,1,,2,]),
                    as.vector(finalscore[,1,,3,]),as.vector(finalscore[,1,,4,]),
                    as.vector(finalscore[,1,,5,]))



# Plot

# Save to png
png(paste0(paths$rst,'/figure1a.png'), width=1000, height=800)
bp=boxplot(two.vectors,outline=FALSE, #col=2:5,
           las = 1,  boxlwd=2, medlwd=2, col=cols,staplewex=0, boxwex=0.8,cex.axis=1.5,cex.names=1.5,#srt = 45,
           names = c('Low','Medium','High','SV Low-Med','SV Low-High')) #,'RW-MRF')
abline(h=1,lwd=1,col='black',lty=3)
dev.off()


# Plot (b) By DGP ----------------------------------------------------------------------------------

# Data
two.vectors = cbind(as.vector(finalscore[,1,,,1]),as.vector(finalscore[,1,,,2]),
                    as.vector(finalscore[,1,,,3]),as.vector(finalscore[,1,,,4])) #,

# Plot
png(paste0(paths$rst,'/figure1b.png'), width=1000, height=800)

bp=boxplot(two.vectors,outline=FALSE, #col=2:5,
           las = 1,  boxlwd=2, medlwd=2, col=cols,staplewex=0, boxwex=0.8,cex.axis=1.5,cex.names=1.5,
           names = c('S1','S2','S3','S4') #,'RW-MRF')
)
abline(h=1,lwd=1,col='black',lty=3)
dev.off()






