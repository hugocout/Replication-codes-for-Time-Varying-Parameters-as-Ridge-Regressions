make_reg_matrix <- function(y,Y,factors,h,ly,lf){
  # Author: Stephane Surprenant
  # Creation: 12/02/2018
  #
  # Description: This function creates a regression matrix
  # containing the dependent variable, y, its lagged values
  # from h to max_y and lagged exogenous regressors, 
  # from lag h to max_f
  #
  # NOTE: y and factors must of same time dimension.
  #
  # OUTPUT
  # First column is dependent variable. All others are
  # regressors. 
  max_y=h+ly-1
  max_f=h+lf-1
  
  bigtY <- nrow(as.matrix(Y))       # Time dimension
  bigtF <- nrow(as.matrix(factors))       # Time dimension
  bign <- ncol(as.matrix(factors)) # Number of factors
  
  lags <- sapply(h:max_y, function(i) c(array(NA,dim=i),Y[1:(bigtY-i)]))
  f=c()
  if(lf>0){
  f <- do.call(cbind, lapply(h:max_f, function(i) 
    rbind(array(NA,dim=c(i,bign)),
          as.matrix(factors[1:(bigtF-i),]))))
  }
  reg <- cbind(y,lags,f)
  return(reg)
}

make_last <- function(y,Y,factors,h,ly,lf){
  # Author: Stephane Surprenant
  # Creation: 13/02/2018
  #
  # Description: This function creates the last values
  # observed corresponding to the regressor matrix made
  # using make_reg_matrix.
  #
  # NOTE: 
  #
  # OUTPUT
  # last is contains last observed values of y and factors
  # in this fashion: y Ly L^2y ..., f1, f2, ..., fK, ... Lf1, ..
  max_y=h+ly-1
  max_f=h+lf-1
  
  y <- as.matrix(y) # Format
  Y <- as.matrix(Y) # Format
  f <- as.matrix(factors) # Format
  n_ly <- max_y-h+1 # Number of y lags
  n_lf <- max_f-h+1 # Number of factors lags
  bign <- dim(f)[2] # Number of factors
  bigt <- dim(y)[1] # time dimension
  
  Ly <- do.call(cbind, lapply(1:n_ly, 
                              function(i) t(as.matrix(tail(Y,i)[1,]))))
  Lf=c()
  if(lf>0){
  Lf <- do.call(cbind, lapply(1:n_lf, 
                              function(i) t(as.matrix(tail(f,i)[1,]))))
  }
  last <- cbind(Ly,Lf)
  return(last)
}
