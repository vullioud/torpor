#' Get Prediction2
#'
#'The function provides the MR estimates and 95% credible interval boundaries
#'at a defined Ta in normothermic and/or torpid stage.
#'@name get_prediction2
#'@aliases get_prediction2
#'@param mod a fitted model from fit_torpor
#'@param Ta a vector of temperatur for which the prediction should be made
#'@return a data frame with predicted values
#'@export
#'@examples
#'\dontrun{
#'data(test_data2)
#'test_mod <- fit_torpor2(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 98,
#'TLC = 28.88,
#'mod = NULL,
#'fitting_options = list(nc = 1))
#'get_prediction2(mod = test_mod, Ta = 1:20)
#'}
get_prediction2 <- function(mod, Ta){


  # retrieved the posterior
  betat <- mod$sims.list$betat
  betac <-mod$sims.list$betac
  inte <-mod$sims.list$inte
  intc <-mod$sims.list$intc
  intr <-mod$sims.list$intr
  Tt <-mod$sims.list$Tt
  tlc <-mod$sims.list$tlc
  Ym <- mod$sims.list$Ym[1]


  X <- length(Ta)
  Ymeant <- rep(NA, X)
  Y975t <- rep(NA, X)
  Y025t <- rep(NA, X)

  Ymeann <- rep(NA, X)
  Y975n <- rep(NA, X)
  Y025n <- rep(NA, X)

  for(i in 1:X) {
    Ymeant[i] <- stats::median(funtorp2(Ta[i], Tt, intr, intc, betat, betac, Ym))
    Y975t[i] <- stats::quantile(funtorp2(Ta[i],Tt, intr, intc, betat, betac, Ym),0.975)
    Y025t[i] <- stats::quantile(funtorp2(Ta[i],Tt, intr, intc, betat, betac, Ym),0.025)
    Ymeann[i] <- stats::median(funnorm2(Ta[i], inte, betat, Ym))
    Y975n[i] <- stats::quantile(funnorm2(Ta[i], inte, betat, Ym),0.975)
    Y025n[i] <- stats::quantile(funnorm2(Ta[i], inte, betat, Ym),0.025)
  }


  #### if mod$sims.list$G
  out1 <- data.frame(Ta = Ta,
                     group = rep("Torp", length(X)),
                     pred =  Ymeant,
                     upr_95 = Y975t,
                     lwr_95 = Y025t)

  out2 <- data.frame(Ta = Ta,
                     group = rep("Norm", length(X)),
                     pred =  Ymeann,
                     upr_95 = Y975n,
                     lwr_95 = Y025n)


  if(length(mod$mean$G[mod$mean$G>1.5]) == 0){
    out <- out2
  } else if (length(mod$mean$G[mod$mean$G<1.5]) == 0){
    out <- out1
  } else {
    out <- rbind(out1, out2)
  }

  return(out)

}

#' funtorp2
#'
#'fit the mod for torpor bats
#'@name funtorp2
#'@param x a temperature
#'@param Tt turning point T
#'@param intr intercept 1
#'@param intc intercept 2
#'@param betat slope 1
#'@param betac slope 2
#'@param Ym mean of Y to back-transform
#'@return a metabolic value

funtorp2 <- function(x, Tt, intr, intc, betat, betac, Ym) {
  out <- (ifelse(x<Tt,intr +betat*x,intc*exp(betac*x)))*Ym
  return(out)
}


#' funnorm2
#'
#'fit the mod for normotermic bats
#'@name funnorm2
#'@param x a temperature
#'@param inte intercept 1
#'@param betat slope 1
#'@param Ym mean of Y to back transform
#'@return a metabolic value
funnorm2 <- function(x, inte, betat, Ym) {
  out <- (inte +betat*x)*Ym
  return(out)
}

