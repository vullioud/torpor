#'get_Q10
#'
#'The function provides the Q10 temperature coefficient,
#'which represents the metabolic rates ratio measured at 10 °C Ta apart.

#'@name get_Q102
#'@param mod a fitted model from fit_torpor
#'@return the Q10 value
#'@export
get_Q102(mod)


get_Q102 <- function(mod){

  betat <- stats::median(mod$sims.list$betat)
  betac <-stats::median(mod$sims.list$betac)
  inte <-stats::median(mod$sims.list$inte)
  intc <-stats::median(mod$sims.list$intc)
  intr <-stats::median(mod$sims.list$intr)
  Tt <-stats::median(mod$sims.list$Tt)
  tlc <-stats::median(mod$sims.list$tlc)
  Ym <- mod$sims.list$Ym[1]



  out <- round((funtorp2(tlc, Tt, intr, intc, betat, betac, Ym)/
    funtorp2(Tt, Tt, intr, intc, betat, betac, Ym))^(10/(tlc-Tt)), digit = 1)

  if(length(mod$mean$G[mod$mean$G>1.5]) == 0){
    warning("no values in torpor")
    out <- NA
  } else {
    out <- out
  }
  return(out)
}
