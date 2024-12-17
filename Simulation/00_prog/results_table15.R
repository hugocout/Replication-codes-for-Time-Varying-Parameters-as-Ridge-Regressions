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
# LOAD RESULTS
# =============================================================================================================

# ShTVP-R and SHTVP-3G models results -------------------------------------------------------------------------

# List all files that match the pattern
mypattern <- 'output_simul5_230817bands_300obs'
file_list <- list.files(paste0(paths$out,"/table15/"), pattern = mypattern)

# Initialize an array to store results
results_array <- array(NA, dim = c(20, 4, 3, 3, 5))

# position in the matrix
modelpos = c(6,7) 
saveall <- list()

for(model in c(1,2)){
  
  model_name <- c("ShTVP-R","ShTVP-3G")[model]
  
  # Processing loop to import data and extract results
  j <- 0
  for (file_name in file_list) {
    j <- ifelse(j > 19, 0, j + 1)  # Reset j after 20 to start over (0-index correction)
    
    load(file.path(paste0(paths$out,"/table15/"), file_name))  # Load the .Rdata file
    
    # Example of processing and capturing specific results, adjust indexing as needed
    selected_scores <- scores[c(1, 10:19, 2, 20, 3:9)[j], seq(modelpos[model], 56, by = 8)[4:7], 1:3, 1:3, 1:5]
    results_array[j, , , , ] <- selected_scores

  }
  
  # Compute means of the results, omitting NA values and rounding
  mean_results <- round(apply(results_array, c(2, 3, 4, 5), mean, na.rm = TRUE), 3)
  saveall[[model]] <- mean_results
}

# 2SRR model results ---------------------------------------------------------------------------------------------

file_name <- "output_simul5_230817bands_2ssr.RData"
load(paste0(paths$out,"/table15/",file_name))
so = apply(scores,c(2,3,4,5),mean)
srrResults <- list()
for(DD in 1:2){
  allo=so[,DD,2,1:5]
  srrResults[[DD]] <- t(c(round(allo*100,1)[26,]/100,round(allo*100,1)[42,]/100,round(allo,2)[50,4:5]))
  
}

# =============================================================================================================
# RESULTS TABLE
# =============================================================================================================

allresults <- cbind(c(srrResults[[1]]),
                    c(saveall[[1]][1, 1, 2, 1:5],saveall[[1]][3, 1, 2, 1:5],saveall[[1]][4, 1, 2, 4:5]),
                    c(saveall[[2]][1, 1, 2, 1:5],saveall[[2]][3, 1, 2, 1:5],saveall[[2]][4, 1, 2, 4:5]),
                    c(srrResults[[2]]),
                    c(saveall[[1]][1, 2, 2, 1:5],saveall[[1]][3, 2, 2, 1:5],saveall[[1]][4, 2, 2, 4:5]),
                    c(saveall[[2]][1, 2, 2, 1:5],saveall[[2]][3, 2, 2, 1:5],saveall[[2]][4, 2, 2, 4:5]))
allresults[1:10,] <- allresults[1:10,]*100
colnames(allresults) <- c("K=6, 2SRR","K=6, ShTVP-R", "K=6, ShTVP-3G", "K=20, 2SRR","K=20, ShTVP-R", "K=20, ShTVP-3G")
rownames(allresults) <- c("68%, low", "med", "high", "SV Low-Med","SV Low-High",
                          "95%, low", "med", "high", "SV Low-Med","SV Low-High",
                          "Cor IQR%, low", "high")
print(allresults)

# Save the results to a .csv file
write.csv(allresults, file = "30_results/table15.csv", row.names = TRUE)
