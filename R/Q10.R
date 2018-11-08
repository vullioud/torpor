#'get_Q10
#'
#'The function provides the Q10 temperature coefficient,
#'which represents the metabolic rates ratio measured at 10 °C Ta apart.

#'@name get_Q10
#'@param mod a fitted model from fit_torpor
#'@return the Q10 value
#'@export
get_Q10 <- function(mod){

  beta1<- stats::median(mod$sims.list$beta1)
  beta2 <- stats::median(mod$sims.list$beta2)
  int1 <- stats::median(mod$sims.list$int1)
  int2 <- stats::median(mod$sims.list$int2)
  int3 <- stats::median(mod$sims.list$int3)
  Tmin <- stats::median(mod$sims.list$Tmin)
  TLC <- 	stats::median(mod$sims.list$tlc)

  out <- round((funtorp(TLC, Tmin, int2, int3, beta1, beta2)/
                  funtorp(Tmin, Tmin, int2, int3, beta1, beta2))
               ^(10/(TLC-Tmin)),digits=1)

  if(length(mod$mean$G[mod$mean$G>1.5]) == 0){
    stop("no values in torpor")
  } else if (length(mod$mean$G[mod$mean$G<1.5]) == 0){
    out <- out
  } else {
    out <- out
  }
return(out)
}
