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
# Table 1
# =============================================================================================================

# This imports all the results data sets of a given folder

mypath = paste0(paths$out,"/table1/")
mypattern = '.Rdata'

tmp.list.1<-list.files(mypath, pattern=mypattern)
tmp.list.2<-list(length=length(tmp.list.1))

# Store results
allresults = array(NA,dim=c(20,4,3,5))

table <- matrix(NA, nrow = 18, ncol = 7)
colnames(table) = c(" ","K=20, 2SRR","K=20, shTVP-R", "K=20, shTVP-3G", "K=50, 2SRR","K=50, shTVP-R", "K=50, shTVP-3G")
ks = c("k/k = 0.2","k/k = 0.5","k/k = 1")
name <- c("k/k = 0.2","Low","Medium","High","SV Low-Med","SV Low-High",
          "k/k = 0.5","Low","Medium","High","SV Low-Med","SV Low-High",
          "k/k = 1","Low","Medium","High","SV Low-Med","SV Low-High")
table[,1] <- name
pos <- c(2:4)

# Make table
for(whichpanel in c(1,2)){
  for (i in 1:length(tmp.list.1)){
    
    load(paste(mypath,tmp.list.1[i],sep=''))
    allresults[i,,,] = scores[c(1,10:19,2,20,3:9)[i],c(2,3,6,7),whichpanel+1,,1:5]
  }
  
  alt <- c()
  for(jiji in 1:3){
    alt <- rbind(alt,rbind(c(NA,NA,NA),t(round(apply(allresults,c(2,3,4),mean,na.rm=TRUE),3)[-2,jiji,])))
  }
  table[1:18,pos] = alt
  pos <- pos + 3
}

# Print results
print(table)

# Save results to csv
write.csv(table, file = paste0(paths$rst,"/table1.csv"), row.names = FALSE)
