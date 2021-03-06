
# Function for clean_datasets using SpATS

library(SpATS)
library(tidyverse)
library(data.table)
source("https://raw.githubusercontent.com/AparicioJohan/SpATS.plus/master/All_additional.R")     # SpATS PLUS

Clean_SpATS <- function(Response, 
                        Geno , 
                        Num_desv=3, 
                        Show_results=TRUE,
                        data=NULL , 
                        fixed=NULL,
                        random=NULL, 
                        col="col", row="row", 
                        name="Exp", iter=NULL){
  
  Datos <- droplevels(data)
  
  w <- 1                                         # Starter 
  k <- Num_desv                                  # Number of standard deviations to consider extreme outliers 
  
  if(inherits(fixed, "character"))
    fixed <- as.formula(fixed)
  if(inherits(random, "character"))
    random <- as.formula(random)
  
  remFinal <- data.frame(Response=as.numeric(),  Genotype = as.character(), col=as.numeric(),row=as.numeric() )
  history <- data.frame(Exp=as.character(), VarE=as.numeric(), VarG=as.numeric(), rep=as.numeric(),r2=as.numeric(), H=as.numeric())
  
  while (w>=1) {
    
    Datos[,Geno] <- as.factor(Datos[,Geno])
    
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
    p <- which( abs(vect_res) > abs( k * sqrt(Var_resi) ) ) # [which.max( abs( vect_res[which(abs(vect_res) > abs(k * sqrt(Var_resi)))] ) )]
    
    remTMP <- data.frame(Response=Datos[p,Response],  Genotype = Datos[p,Geno], col=Datos[p,"col"],row=Datos[p,"row"] )
    remFinal <- rbind(remTMP,remFinal)
    
    Datos[p,Response] <- NA
    
    # History
    VarG <- as.numeric(Modelo$var.comp[Geno])
    VarE <- as.numeric(Modelo$psi[1])
    replicate <- NA
    r2 <- as.numeric(R.square(Modelo))
    Sum <- data.frame(Exp= paste0(name) , VarE=VarE, VarG=VarG, rep=replicate , r2=r2)
    Sum$H <- as.numeric(getHeritability(Modelo))      
    history <- rbind(history,Sum)               
    
    if (Show_results==TRUE) {
      cat("\n" , "N_xtreme_residuals:" , w, "\tHeritability:", getHeritability(Modelo), "\tAIC:", AIC(Modelo)  )  
    }
    
    if(!is.null(iter)) w <- 0
                   
  }
  
  Clean_VEF <- Datos
  # Clean_VEF <- Clean_VEF[,-c(ncol(Clean_VEF),ncol(Clean_VEF)-1 )]
  
  BLUPs <- predict(object = Modelo, which = Geno) %>% 
    dplyr::select(all_of(Geno), predicted.values, standard.errors) %>%
    mutate(predicted.values=predicted.values-mean(predicted.values,na.rm=T))
                   
  names(BLUPs) <- c("level","estimate","std.error") 
  
  # names(BLUPs)[ncol(BLUPs)] <- "Trial"
  VarG <- as.numeric(Modelo$var.comp[Geno])
  VarE <- as.numeric(Modelo$psi[1])
  replicate <- NA
  r2 <- as.numeric(R.square(Modelo))
  
  Sum <- data.frame(Exp= paste0(name) , VarE=VarE, VarG=VarG, rep=replicate , r2=r2)
  Sum$H <- as.numeric(getHeritability(Modelo))
  
  k=list(BLUPs=BLUPs, data_clean=Clean_VEF, Model=Modelo, Remove=remFinal, Summ=Sum, History=history)
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
