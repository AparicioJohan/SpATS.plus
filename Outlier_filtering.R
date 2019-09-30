
# Function for clean_datasets using SpATS

library(SpATS)
library(tidyverse)
library(data.table)
source("https://raw.githubusercontent.com/AparicioJohan/SpATS.plus/master/All_additional.R")     # SpATS PLUS

Clean_SpATS <- function(Response, Geno , Num_desv=3, Show_results=TRUE, data=NULL , fixed=NULL, random=NULL, col="col", row="row"){
  
  Datos <- droplevels(data)
  
  w <- 1                                         # Starter 
  k <- Num_desv                                  # Number of standard deviations to consider extreme outliers 
  
  if(inherits(fixed, "character"))
    fixed <- as.formula(fixed)
  if(inherits(random, "character"))
    random <- as.formula(random)
  
  while (w>=1) {
    
    Datos[,Geno]<- as.factor(Datos[,Geno])
    
    Datos$col_f = factor(Datos[,col])
    Datos$row_f = factor(Datos[,row])
    
    if(!all(sapply(all.vars(random), function(x, data) is.factor(data[,x]), data = Datos))) {
      stop("All variables indicated in argument 'random' should be factors")
    }
    
    ncols = length(unique(Datos[,col]))
    nrows = length(unique(Datos[,row]))
    
    Modelo = SpATS(response=Response,
                   genotype=Geno, genotype.as.random=T,
                   fixed=fixed ,
                   spatial = ~ PSANOVA(col, row, nseg = c(ncols,nrows), degree=c(3,3),nest.div=2),
                   random =random , data=Datos,
                   control = list(tolerance=1e-03, monitoring=0))
    
    Obs <- Modelo$nobs
    Eff.dim <- sum(c(Modelo$eff.dim))
    Var_resi <- sum( na.omit(residuals(Modelo))^2 ) / (Obs - Eff.dim) # Sum(Errores^ 2)/ED_e
    vect_res <- residuals(Modelo)
    
    # Number of extreme residuals (above k standard deviations) in this iteration
    w <- length( which( abs(vect_res) > abs(k * sqrt(Var_resi)) ) )
    
    # What is the most extreme residual ?
    p <- which( abs(vect_res) > abs(k * sqrt(Var_resi)) )[which.max( abs( vect_res[which(abs(vect_res) > abs(k * sqrt(Var_resi)))] ) )]
    
    Datos[p,Response] <- NA
    
    
    if (Show_results==TRUE) {
      cat("\n" , "N_xtreme_residuals:" , w, "\tHeritability:", getHeritability(Modelo), "\tAIC:", AIC(Modelo), "\n"  )  
    }
    
    
  }
  
  Clean_VEF <- Datos
  Clean_VEF <- Clean_VEF[,-c(ncol(Clean_VEF),ncol(Clean_VEF)-1 )]
  
  BLUPs <- predict(object = Modelo, which = Geno) %>% 
      select(Geno, predicted.values)

  
  names(BLUPs)[ncol(BLUPs)] <- "Trial"
  cat('\n')
  
  k=list(BLUPs=BLUPs, data_clean=Clean_VEF)
  k
  
}

#################################
##          Example 
#################################


# temp <-  Clean_SpATS(Response = "yield", 
#              Geno = "geno", 
#              Num_desv = 3 , Show_results = T,
#              data = wheatdata, 
#              fixed = NULL,
#              random = ~ col_f + row_f,
#              col= "col" , row = "row"  )