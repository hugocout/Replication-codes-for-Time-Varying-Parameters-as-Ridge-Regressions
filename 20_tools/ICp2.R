ICp2 <- function(mse, n_fac, bigt, bign){
  
  # Author: Stephane Surprenant
  # Date: 13/10/2017
  
  # This function evaluates the Bai & Ng (2002) information criteria for
  # choosing the number of factors in a factor model using the mse of a
  # previously computed factor model and the corresponding number of factors,
  # n_fac. The computation also requires the temploral and cross-sectional
  # dimensions, bigt and bign.
  #
  # mse = V(k,\hat{F}^k) = 1/NT sum_i sum_t \hat{e}_t^2 
  # C_{NT} = min{sqrt(N), cqrt(T)}
  # IC_{p2} = ln(V(k,\hat{F}^k)) + k\frac{N+T}{N} ln C_{NT}^2
  #
  # Bai, J. and Serena Ng. 2002. Determining the number of factors in
  #     approximate factor models. Econometrica, 70(1), p. 191-221.
  #
  # ========================================================================= #
  #                               ESTIMATION
  # ========================================================================= #
  
  k <- n_fac
  v <- mse
  
  c_nt <- min(c(sqrt(bign), sqrt(bigt)))
  CT <- k*((bign+bigt)/(bign*bigt))*log(c_nt^2) # Penalty function
  
  icp2 <- log(v) + CT
  
  return(icp2)
}