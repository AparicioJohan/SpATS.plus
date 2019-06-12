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