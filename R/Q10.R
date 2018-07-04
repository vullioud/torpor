#'get_Q10
#'
#'compute the Q10 values
#'@name get_Q10
#'@param model a fitted model from fit_torpor
#'@return the Q10 value

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
return(out)
}
