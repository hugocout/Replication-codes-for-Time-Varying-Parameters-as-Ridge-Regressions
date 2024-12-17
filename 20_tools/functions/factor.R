factor <- function(data,n_fac,id=1){
  
  # Author: Stephane Surprenant
  # Creation: 12/10/2017
  #
  # This function estimates factors and loadings using principal component
  # analysis. The model is X = F^0 Lambda^0' + e, with X is T by N data matrix,
  # F is the the T by r factor matrix, Lambda is the N by r loadings matrix and
  # e is the T by r errors matrix. (Supercript k indicates k factors used.)
  #
  # Under F^k'F^k/T = I, maximizing the trace of F^k'(XX')F^k, \hat{F}^k = 
  # sqrt(T) times the eigenvectors of the k largest eigenvalues of XX'. Under 
  # Lambda'Lambda/N = I, we have Lambda_k_hat = sqrt(N) times the eignevectors
  # of X'X. It also implies \hat{F}^k = X \hat{Lambda}^k/N. (Bai & Ng, 2002).
  # The second choice is implemented here.
  #
  # Bai, J. and Serena Ng. 2002. Determining the number of factors in
  #     approximate factor models. Econometrica, 70(1), p. 191-221.
  #
  # ========================================================================= #
  #                               ESTIMATION
  # ========================================================================= #
  
  X <- data  # Reassign
  r <- n_fac # Reassign
  
  bign <- dim(X)[2]
  bigt <- dim(X)[1]
  
  svd <- svd(crossprod(X))         # svd$u = svd$v (symmetry); svd$d = eigenvalues
  lambda <- svd$u[,1:r]*sqrt(bign) # r th column times r th biggest eigenvalue
  f_hat <- X%*%lambda/bign         # factors
  e_hat <- X - f_hat%*%t(lambda)   # errors
  mse <- sum(e_hat^2)/(bign*bigt)
  
  if(id!=3){lambda = lambda %*% diag(r)} #does nothing
  
  if(id==3){
    F.F = pracma::sqrtm(crossprod(f_hat))

    lambda =  lambda %*% F.F$B
    f_hat = f_hat%*%F.F$Binv

    #rotation to PC3 identification
    lambda=lambda %*% solve(lambda[1:r,])
    f_hat = tcrossprod(f_hat,lambda[1:r,])
  }
  
  results <- list(factors = f_hat, lambda = lambda, mse = mse)
  
  return(results)
}
