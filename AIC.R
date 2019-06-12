AIC.SpATS <- function(Modelo){
  p <- length(Modelo$coeff[!attr(Modelo$coeff,"random")])   # number of fixed coefficients
  q <- length(Modelo$var.comp) + 1                          # number of random components + 1 (residual variance)
  lo.lik <- Modelo$deviance/-2                              # log Likelihood 
  aic <- -2*lo.lik + 2*(p+q)                                # Akaike Information Criterion
  return(aic)
}
