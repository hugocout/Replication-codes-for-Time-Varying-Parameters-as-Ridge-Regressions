
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

# =============================================================================================================
# SIMULATION PARAMETERS
# =============================================================================================================

# Load the simulation types
Stype = 5#all_options$H[it_pos] 
obs <- 301

# Parameters
whichsimul = paste('simul',Stype,sep='')
endofsample = 0*15
setups = c(2,1,3,4,5)
crit = 2
models = 4
totcol = 15
fracvec = c(.8,.5,0) 
DDvec = c(6,6,6)#c(6,20,100)

# =============================================================================================================
# SIMULATION
# =============================================================================================================

for(DD in c(2)){
  for(frac in c(3)){
    for(setup in 4){
      pos <- 1
      betas <- matrix(NA, nrow = 4, ncol = obs-1)
      for(Stype in c(0,3,4,5)) {
        whichsimul = paste('simul',Stype,sep='')
        
        T=obs
        L<-1
        Mbase<-DDvec[DD]
        compressor = .5
        
        # Setups ------------------------------------------------------
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
        s=100
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
        
        # Plot all betas
        betas[pos,] <- beta[1,]
        
        pos <- pos+1
      }
      
      png(paste0(paths$rst,'/figure4.png'), width = 1000, height = 800)
      # Plot all betas
      title <- paste0('DD',DD,', frac',frac,', setup',setup)
      plot(betas[1,], lwd = 2, type = 'l', ylim = c(min(betas),max(betas)), main = title)
      lines(betas[2,], lwd = 2, col = 'red')
      lines(betas[3,], lwd = 2, col = 'blue')
      lines(betas[4,], lwd = 2, col = 'green')
      legend('topleft', legend = paste0("f",1:4), col = c('black', 'red', 'blue', 'green'), lwd = 2,
             cex = 1.5)
      dev.off()
    }
  }
}