#' get Prediction
#'
#' The function [get_prediction()] provides the predicted MR and 95CI bounds at
#' a given Ta, in normothermic and/or torpid stage.
#'
#'@name get_prediction
#'@aliases get_prediction
#'@param mod a fitted model from the function [fit_torpor()]
#'@param Ta a vector of temperatur for which the prediction should be made
#'@return a data frame with predicted values
#'@export
#'@examples
#'\dontrun{
#'data(test_data2)
#'test_mod <- fit_torpor(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 1.49,
#'TLC = 28.88,
#'fitting_options = list(nc = 1))
#'get_prediction(mod = test_mod, Ta = 1:20)
#'}
get_prediction <- function(mod, Ta){


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
    Ymeant[i] <- stats::median(funtorp(Ta[i], Tt, intr, intc, betat, betac, Ym))
    Y975t[i] <- stats::quantile(funtorp(Ta[i],Tt, intr, intc, betat, betac, Ym),0.975)
    Y025t[i] <- stats::quantile(funtorp(Ta[i],Tt, intr, intc, betat, betac, Ym),0.025)
    Ymeann[i] <- stats::median(funnorm(Ta[i], inte, betat, Ym))
    Y975n[i] <- stats::quantile(funnorm(Ta[i], inte, betat, Ym),0.975)
    Y025n[i] <- stats::quantile(funnorm(Ta[i], inte, betat, Ym),0.025)
  }


  #### if mod$sims.list$G
  out1 <- data.frame(Ta = Ta,
                     group = rep("Torpor", length(X)),
                     pred =  Ymeant,
                     upr_95 = Y975t,
                     lwr_95 = Y025t)

  out2 <- data.frame(Ta = Ta,
                     group = rep("Euthermy", length(X)),
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

#' funtorp
#'
#'fit the mod for torpor bats // used internally only
#'@name funtorp
#'@param x a temperature
#'@param Tt turning point T
#'@param intr intercept 1
#'@param intc intercept 2
#'@param betat slope 1
#'@param betac slope 2
#'@param Ym mean of Y to back-transform
#'@return a metabolic value

funtorp <- function(x, Tt, intr, intc, betat, betac, Ym) {
  out <- (ifelse(x<Tt,intr +betat*x,intc*exp(betac*x)))*Ym
  return(out)
}


#' funnorm
#'
#'fit the mod for normotermic bats // used internally only
#'@name funnorm
#'@param x a temperature
#'@param inte intercept 1
#'@param betat slope 1
#'@param Ym mean of Y to back transform
#'@return a metabolic value
funnorm <- function(x, inte, betat, Ym) {
  out <- (inte +betat*x)*Ym
  return(out)
}

################################################################################
#' Get classification
#'
#'The function get_classification() returns the raw data with the related
#'predicted stage values (between 1 and 2), which leads to the stage
#'classification (torpor or euthermy). Additionally, it also provides the
#'predicted MR at the given Ta.
#'
#'@name get_classification
#'@aliases get_classification
#'@param mod a fitted model from fit_torpor
#'@return a data frame with classification and predicted MR
#'@export
#'@examples
#'\dontrun{
#'data(test_data2)
#'test_mod <- fit_torpor(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 1.49,
#'TLC = 28.88,
#'fitting_options = list(nc = 1))
#'get_classification(mod = test_mod)
#'}

get_classification <- function(mod){
x <- data.frame(measured_MR = (mod$data$Y)*mod$data$Ym[1],
                measured_Ta = mod$data$Ta,
                predicted_stage = mod$mean$G)

betat <- mod$sims.list$betat
betac <-mod$sims.list$betac
inte <-mod$sims.list$inte
intc <-mod$sims.list$intc
intr <-mod$sims.list$intr
Tt <-mod$sims.list$Tt
tlc <-mod$sims.list$tlc
Ym <- mod$sims.list$Ym[1]


X <- nrow(x)
x$predicted_MR <- rep(NA, X)
x$classification <- ifelse(x$predicted_stage > 1.5, "Torpor", "Euthermy")

for(i in 1:nrow(x)) {
  if (x$predicted_stage[i] > 1.5) {
  x$predicted_MR[i] <- stats::median(funtorp(x$measured_Ta[i], Tt, intr, intc, betat, betac, Ym))
  } else {
  x$predicted_MR[i] <- stats::median(funnorm(x$measured_Ta[i], inte, betat, Ym))
  }
}
return(x)
}
