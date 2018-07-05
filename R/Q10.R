#'get_Q10
#'
#'compute the Q10 values
#'@name get_Q10
#'@param model a fitted model from fit_torpor
#'@return the Q10 value
#'@export

get_Q10 <- function(model){

  beta1<- stats::median(model$sims.list$beta1)
  beta2 <- stats::median(model$sims.list$beta2)
  int1 <- stats::median(model$sims.list$int1)
  int2 <- stats::median(model$sims.list$int2)
  int3 <- stats::median(model$sims.list$int3)
  Tmin <- stats::median(model$sims.list$Tmin)
  TLC <- 	stats::median(model$sims.list$tlc)

  out <- round((funtorp(TLC, Tmin, int2, int3, beta1, beta2)/
                  funtorp(Tmin, Tmin, int2, int3, beta1, beta2))
               ^(10/(TLC-Tmin)),digits=1)

  if(length(model$mean$G[model$mean$G>1.5]) == 0){
    warning("error no values in torpor")
    out <- NA
  } else if (length(model$mean$G[model$mean$G<1.5]) == 0){
    out <- out
  } else {
    out <- out
  }
return(out)
}
