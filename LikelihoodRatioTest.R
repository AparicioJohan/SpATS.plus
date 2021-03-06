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
