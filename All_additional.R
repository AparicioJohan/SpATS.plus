###########################
# ALL additional function
###########################

library(SpATS)


# AIC 
AIC.SpATS <- function(Modelo){
  p <- length(Modelo$coeff[!attr(Modelo$coeff,"random")])   # number of fixed coefficients
  q <- length(Modelo$var.comp) + 1                          # number of random components + residual variance
  lo.lik <- Modelo$deviance/-2                              # log Likelihood 
  aic <- -2*lo.lik + 2*(p+q)                                # Akaike Information Criterion
  return(aic)
}

# BIC
BIC.SpATS <- function(Modelo){
  p <- length(Modelo$coeff[!attr(Modelo$coeff,"random")])   # number of fixed coefficients
  q <- length(Modelo$var.comp) + 1                          # number of random components + residual variance
  lo.lik <- Modelo$deviance/-2                              # log Likelihood 
  n <- Modelo$nobs                                          # number of observations
  bic <- -2*lo.lik + log(n)*(p+q)                           # Bayesian information criterion
  return(bic)
}

# Likelihood Ratio Test
Lik.ratio.test <- function(Model_nested,Model_full){
  
  p1 <- length(Model_nested$coeff[!attr(Model_nested$coeff,"random")])
  p2 <- length(Model_full$coeff[!attr(Model_full$coeff,"random")])
  if(p1!=p2) stop("\nFixed factors should be the same!\n")
  
  lo.lik1 <- Model_nested$deviance/-2
  lo.lik2 <- Model_full$deviance/-2
  
  d <- 2*(lo.lik2-lo.lik1)
  
  r1 <- length(Model_nested$var.comp) + 1
  r2 <- length(Model_full$var.comp) + 1
  if(r1==r2) stop("\nEqual number of variance components")
  
  p.value1=round(1 - pchisq(d, r2-r1), 3)
  
  siglevel <- 0
  if (abs(p.value1) < 0.05) {
    siglevel <- "*"
  }
  else {
    siglevel <- "Not signif"
  }
  if (abs(p.value1) < 0.01) {
    siglevel <- "**"
  }
  if (abs(p.value1) < 0.001) {
    siglevel <- "***"
  }
  
  names(p.value1) <- "p-value"
  cat("\n========================================================")
  cat("\n", "Model 1 ==> ", "AIC:", AIC(Model_nested),"|",  "BIC:", BIC(Model_nested),"|", "Varcomp:" ,r1)
  cat("\n", "Model 2 ==> ", "AIC:", AIC(Model_full),"|", "BIC:", BIC(Model_full),"|", "Varcomp:" ,r2, "\n")
  cat("-----------------------------\n")
  cat("Lower AIC and BIC is better model.\n\n")
  cat("\n", "Likelihood Ratio Test")
  cat("\n", "p-value =" ,p.value1, siglevel)
  cat("\n", "Sig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n")
  cat("========================================================")
  
}

# Pseudo R square
R.square <- function(Model){
  
  response      <- Model$data[,Model$model$response]  
  mean.response <- mean(response,na.rm = T)
  fitted        <- Model$fitted
  SS.fitted     <- sum( (response-fitted)^2 , na.rm = T)
  SS.response   <- sum( (response-mean.response)^2 , na.rm = T )
  
  R <- 1- SS.fitted/SS.response
  
  names(R) <- "r.square"
  return(round(R,3))
}


###########################
#       Example
###########################


# data(wheatdata)
# summary(wheatdata)
# 
# # Create factor variable for row and columns
# wheatdata$R <- as.factor(wheatdata$row)
# wheatdata$C <- as.factor(wheatdata$col)
# 
# m0 <- SpATS(response = "yield", spatial = ~ PSANOVA(col, row, nseg = c(10,20), degree = 3, pord = 2), 
#             genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
#             control =  list(tolerance = 1e-03))
# 
# m1 <- SpATS(response = "yield", spatial = ~ PSANOVA(col, row, nseg = c(10,20), degree = 3, pord = 2), 
#             genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C + rep, data = wheatdata, 
#             control =  list(tolerance = 1e-03))
# 
# Lik.ratio.test(Model_nested = m0,Model_full = m1)
# 
# AIC(m0); AIC(m1)
# 
# BIC(m0); BIC(m1)
# 
# R.square(m0);R.square(m1)


