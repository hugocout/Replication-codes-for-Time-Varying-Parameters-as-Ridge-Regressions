EM_sw <- function(data, n, it_max=50){
  #
  # Author: Stephane Surprenant
  # Creation: 16/11/2017
  #
  # Description:
  # This function inputs missing values in 'data' by using an
  # Expectation-Maximization algorithm. The data must is standerdized and
  # centered at zero. The algorithm begins by replacing the missing values
  # by 0 (the unconditional mean of the available observations). Then, it
  # estimates a factor model through principal component. Missing values
  # are then replaced by the predicted common component of the series.
  #
  # INPUTS: data: a matrix of data with time series as column vectors;
  #         n: number of factors considered
  # OUTPUT: X: a balanced panel
  #
  # Note: Requires package pracma (for repmat) and factor.R for PCA.
  #
  # ========================================================================= #
  #                               ESTIMATION
  # ========================================================================= #
  
  X <- as.matrix(data)
  
  # 1. Standerdize data
  size <- dim(X)
  a <- standard(X)
  X0 <- a$Y
  std <- a$std
  mu <- a$mean
  
  # 2. First replacement and PCA
  X0[is.na(X)] <- 0 # Unconditional mean of standard data is zero
  X0[is.na(X0)] <- 0 # Unconditional mean of standard data is zero
  
  a <- factor(X0, n_fac=n) # PCA
  f_1 <- a$factors
  L_1 <- a$lambda
  x_hat <- X0
  cc <- f_1%*%t(L_1) # Common component 
  x_hat[is.na(X)] <- cc[is.na(X)]
  x_hat <- x_hat*repmat(std,size[1],1) + repmat(mu,size[1],1)
  
  # 3. Initialize algorithm
  it <- 0      # Initialize iterations
  err <- 999   # Initialize error criterion to arbitrarily high value
  # it_max <- 50 # Selection maximum iterations
  con <- 1e-6  # Set convergence criterion
  old_f <- array(0, dim=c(size[1],1))
  
  while(it < it_max && err > con){
    # 1. Standerdize dataset
    a <- standard(x_hat) 
    std_new <- a$std
    mu_new <- a$mean
    x_new <- a$Y
    
    old_f <- f_1 # Save factors from previous iteration
    
    # 2. Extract factor through PCA
    a <- factor(x_new, n_fac=n)
    f_1 <- a$factors
    L_1 <- a$lambda
    mse <- a$mse
    cc <- f_1%*%t(L_1) # Common component 
    x_hat <- x_new
    x_hat[is.na(X)] <- cc[is.na(X)]
    x_hat <- x_hat*repmat(std_new,size[1],1) + repmat(mu_new,size[1],1)
    
    mean_old_f2 <- apply(old_f^2, c(2), mean)
    mean_new_f2 <- apply(f_1^2, c(2), mean)
    
    err <- abs(mean(mean_old_f2 - mean_new_f2))
    
    it <- it + 1;
  }
  x_hat[!is.na(X)] <- X[!is.na(X)]
  results <- list(data=x_hat, factors=f_1, lambda=L_1, iterations=it, mse=mse)
  return(results)
}

# Subfunction to standerdize dataset ======================================== #
standard <- function(Y){
  Y < as.matrix(Y)
  size <- dim(Y)
  
  mean_y <-apply(Y, c(2), mean, na.rm=TRUE)
  sd_y <- apply(Y, c(2), sd, na.rm=TRUE)
  
  Y0 <- (Y - repmat(mean_y, size[1],1))/repmat(sd_y, size[1],1)
  return(list(Y=Y0, mean=mean_y, std=sd_y))
}