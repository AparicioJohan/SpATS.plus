
require(ggplot2,ggpubr)
autoplot.SpATS <- function(Model,...){

  response      <- Model$data[,Model$model$response]          # Response Variable
  genotype      <- Model$data[,Model$model$geno$genotype]     # Genotypes
  var.res       <- Model$psi[1]                               # Residual Variance
  fitted        <- Model$fitted                               # Fitted Values
  residuals     <- Model$residuals                            # Residuals
  stand.res     <- residuals/sqrt(var.res)                    # Scaled Residuals
  sqrt.stand    <- sqrt(abs(stand.res))                       # Scale Location
  inf.s         <- data.frame(genotype,response,fitted,residuals,stand.res,sqrt.stand) 
  
  p1 <- ggplot(inf.s, aes(x=fitted,y=stand.res))+geom_point(...)+ labs(x="Fitted values",y="Standardized residuals", title = "Residuals vs Fitted")
  p2 <- ggplot(inf.s, aes(sample=stand.res))+stat_qq(...)+stat_qq_line()+labs(x="Theorical Quantile",y="Standardized residuals",title = "Normal Q-Q")
  p3 <- ggplot(inf.s, aes(x=fitted,y=sqrt.stand))+geom_point(...)+ labs(x="Fitted values",y=expression(sqrt("|Standardized residuals|")), title = "Scale-Location")
  p4 <- ggplot(inf.s, aes(x=as.numeric(genotype),y=stand.res))+ geom_point(...)+labs(x="Factor Genotype",y="Standardized residuals",title = "Residuals vs Factor levels")
  
  return(suppressWarnings(ggarrange(p1,p2,p3,p4)))
}



