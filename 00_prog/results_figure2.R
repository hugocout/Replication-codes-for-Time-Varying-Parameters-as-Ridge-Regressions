
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
myPKGs <- c('pracma','Hmisc','forecast','RColorBrewer','ggplot2','reshape2','extrafont',
            'grid','ggthemes')

# check to see if each package is installed, and if not add it to a list to install
InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]

# install any needed packages
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

# Libraries
library(pracma)     # For repmat 
library(Hmisc)
library(forecast)   # Diebold-Mariano test
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(extrafont)
library(grid)
library(ggthemes)

# Plot custom functions
theme_Publication <- function(base_size=20, base_family="Arial") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.3), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = 'black'),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(angle=45,vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.6, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.margin=margin(-20,5,5,-20),
            legend.box.margin=margin(-5,-5,-5,-5),
            #legend.text=element_text(size=12),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale(name = "","fill","Publication",
                 manual_pal(values = c("#386cb0",brewer.pal('Dark2',n=3)[1],
                                       "#fdb462","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",
                 manual_pal(values = c("#386cb0","#fdb462",brewer.pal('Dark2',n=3)[1],
                                       "#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

# =============================================================================================================
# LOAD RESULTS
# =============================================================================================================

# Load and create storage
load(paste(paths$out,'/table16to18/TVPfcst_111_bayes.Rdata',sep=''))
dimin=dim(forecast)
dimin[4]=5
old_err = array(NA,dim=c(dimin,4))
benchmarks = array(NA,dim=c(dimin[1],dimin[2],dimin[3],4))

for(mod in 1:2){
  for(v in 1:5){
    for(h in c(1,2,4)){
      
      # Load
      load(paste(paths$out,'/table16to18/TVPfcst_',v,h,mod,'.Rdata',sep=''))
      h_pos <- h
      if(h == 4) {
        h_pos <- dim(forecast)[2]
      }
      
      old_err[,h,v,1:3,mod]=forecast[,h_pos,v,1:3]
      old_err[,h,v,4,mod]=(0.5*forecast[,h_pos,v,2])+(0.5*forecast[,h_pos,v,1])
      
      if(mod == 1) {
        benchmarks[,h,v,1]=forecast[,h_pos,v,1]
        benchmarks[,h,v,2]=forecast[,h_pos,v,1]
      }
      
      # Load block CV results
      load(paste(paths$out,'/table16to18/TVPfcst_',v,h,mod,'_block.Rdata',sep=''))
      if(h == 4) {
        h_pos <- dim(forecast)[2]
      }
      
      old_err[,h,v,5,mod]=forecast[,h_pos,v,2]
      
      if(mod == 1) {
        benchmarks[,h,v,3]=forecast[,h_pos,v,1]
      }
      
      # Load bayes results
      load(paste(paths$out,'/table16to18/TVPfcst_',v,h,mod,'_bayes.Rdata',sep=''))
      if(h == 4) {
        h_pos <- dim(forecast)[2]
      }
      
      old_err[,h,v,3,mod]=forecast[,h_pos,v,2]
      
      if(mod == 1) {
        benchmarks[,h,v,4]=forecast[,h_pos,v,1]
      }
      
      
    }
  }
}
old_err=old_err[158:222,-3,,,] # Remove empty columns
benchmarks=benchmarks[158:222,-3,,] # Remove empty columns

# =============================================================================================================
# RMSES BARPLOTS
# =============================================================================================================

# Compute RMSE --------------------------------------------------------------
mspe_mattemp=sqrt(apply(old_err^2,c(2,3,4,5),mean))
mspe_mat = mspe_mattemp
mspe_benchmarks = sqrt(apply(benchmarks^2,c(2,3,4),mean))
for(m in 1:5){
  for(mod in 1:2){
    
    if(m == 5) {
      mspe_mat[,,m,mod] = mspe_mattemp[,,m,mod]/mspe_benchmarks[,,3] 
    }else if(m == 3) {
      mspe_mat[,,m,mod] = mspe_mattemp[,,m,mod]/mspe_benchmarks[,,4]
    }else{
      mspe_mat[,,m,mod] = mspe_mattemp[,,m,mod]/mspe_benchmarks[,,1] 
    }
  }
}

# Barplots -----------------------------------------------------------------
targets_names <- c("GDP","UR","INF","IR","SPREAD")
horizons <- c(1,2,4)
for(v in c(1:5)){
  for(h in c(1,2,3)){
    
    # Data preparation -------------------------------------------------------
    err=old_err[,h,v,c(1,2,3,4,5),1:2]
    alpha=0.1
    allthings = mspe_mat[h,v,,1:2]
    dimnames(allthings)[[2]] = c('AR','ARDI')
    dimnames(allthings)[[1]] = c('Plain','2SRR','Bayes','2SRR \n H&H','2SRR \n BCV')
    tp.long = reshape2::melt(allthings[,])
    colnames(tp.long)=c('model','rsq') 
    
    
    
    # Diebold-Mariano test ---------------------------------------------------
    allthings = array(NA, dim =c(5,2))
    for(m in 1:5){
      if(m == 5) {
        err_base=benchmarks[,h,v,3]
      }else if(m == 3){
        err_base=benchmarks[,h,v,4]
      }else{
        err_base=benchmarks[,h,v,1]
      }
      for(mod in 1:2){
        if(m==1 & mod==1){
          allthings[m,mod]=0
        }else{
          dm_test <- base::ifelse(I(dm.test(err[,m,mod],err_base,h=c(1,2,4)[h],alternative='two.sided',power=2)$p.value<alpha),0,-1)
          allthings[m,mod]=dm_test+2
        }
      }
    }
    dimnames(allthings)[[2]] = c('AR','ARDI')
    dimnames(allthings)[[1]] = c('Plain','2SRR','Bayes','2SRR \n H&H','2SRR \n BCV')
    
    # Make plot ---------------------------------------------------------------
    tp.long2 = reshape2::melt(allthings[,])
    colnames(tp.long2)=c('model','test')
    
    data_final = as.data.frame(cbind(tp.long,tp.long2[,3]))
    colnames(data_final)=c('model','cat','rsq','test') 
    data_final=data_final[complete.cases(data_final),]
    
    plot_title = paste(targets_names[v], ', h=',horizons[h],sep='')
    p=ggplot(data = data_final, aes(x = model, y = rsq, fill=as.factor(test))) + 
      geom_bar(stat='identity',color='black') + 
      theme_Publication()+
      facet_grid(~cat,scales='free',space='free_x')+
      xlab('')+ylab('')+labs(color='')+
      theme(strip.text = element_text(face="bold", colour = "white",size=25,family="Arial"), #
            strip.background=element_rect(colour="black",fill="black"))+
      geom_hline(data = data.frame(yint=data_final[data_final$model=='Plain'&data_final$cat=='AR',3],cat="AR"), aes(yintercept = yint), linetype = "dashed",lwd=0.75) + 
      geom_hline(data = data.frame(yint=data_final[data_final$model=='Plain'&data_final$cat=='ARDI',3],cat="ARDI"), aes(yintercept = yint), linetype = "dashed",lwd=0.75) +
      theme(legend.position = "none")+scale_fill_Publication() +
      scale_fill_manual(breaks = c(0, 1, 2),
                        values=c("#386cb0","#fdb462",brewer.pal('Dark2',n=3)[1]))+
      ggtitle(plot_title)
    
    # Print
    plot(p)
  }
}
