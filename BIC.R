BIC.SpATS <- function(Modelo){
  p <- length(Modelo$coeff[!attr(Modelo$coeff,"random")])   # number of fixed coefficients
  q <- length(Modelo$var.comp) + 1                          # number of random components + residual variance
  lo.lik <- Modelo$deviance/-2                              # log Likelihood 
  n <- Modelo$nobs                                          # number of observations
  bic <- -2*lo.lik + log(n)*(p+q)                           # Bayesian information criterion
  return(bic)
}