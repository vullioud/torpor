#' Get Prediction2
#'
#'The function provides the MR estimates and 95% credible interval boundaries
#'at a defined Ta in normothermic and/or torpid stage.
#'@name get_prediction2
#'@aliases get_prediction2
#'@param model a fitted model from fit_torpor
#'@param Ta a vector of temperatur for which the prediction should be made
#'@return a data frame with predicted values
#'@examples
#'\dontrun{
#'data(test_data2)
#'test_mod <- fit_torpor2(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 98,
#'TLC = 28.88,
#'Model = NULL,
#'fitting_options = list(nc = 1))
#'get_prediction(test_mod, 20)
#'}
#'@export
get_prediction2 <- function(model, Ta){
  X <- Ta


  betat<- model$sims.list$betat
  betac <-model$sims.list$betac
  inte <-model$sims.list$inte
  intc <-model$sims.list$intc
  intr <-model$sims.list$intr
  Tt <-model$sims.list$Tt
  tlc <-model$sims.list$tlc

  ##### tau to check if correct
  model$sims.list$tauy
  sd(model$sims.list$tauy[, 1])
  sd(model$sims.list$tauy[, 2])

  X <- Ta
  Ymeant <- rep(NA, X)
  Y975t <- rep(NA, X)
  Y025t <- rep(NA, X)
  Ymeann <- rep(NA, X)
  Y975n <- rep(NA, X)
  Y025n <- rep(NA, X)

  for(i in 1:length(Ta)) {
    Ymeant[i] <- stats::median(funtorp(X[i], Tt, intr, intc, betat, betac, Ym))
    Y975t[i] <- stats::quantile(funtorp(X[i],Tt, intr, intc, betat, betac, Ym),0.975)
    Y025t[i] <- stats::quantile(funtorp(X[i],Tt, intr, intc, betat, betac, Ym),0.025)
    Ymeann[i] <- stats::median(funnorm(X[i], int1, beta1))
    Y975n[i] <- stats::quantile(funnorm(X[i], int1, beta1),0.975)
    Y025n[i] <- stats::quantile(funnorm(X[i], int1, beta1),0.025)
  }


  #### if model$sims.list$G
  out1 <- data.frame(Ta = X,
                     group = rep("Torp", length(X)),
                     pred =  Ymeant,
                     upr_95 = Y975t,
                     lwr_95 = Y025t)

  out2 <- data.frame(Ta = X,
                     group = rep("Norm", length(X)),
                     pred =  Ymeann,
                     upr_95 = Y975n,
                     lwr_95 = Y025n)


  if(length(model$mean$G[model$mean$G>1.5]) == 0){
    out <- out2
  } else if (length(model$mean$G[model$mean$G<1.5]) == 0){
    out <- out1
  } else {
    out <- rbind(out1, out2)
  }

  return(out)

}

#' funtorp
#'
#'fit the model for torpor bats
#'@name funtorp
#'@param x a temperature
#'@param Tmin turning point T
#'@param int2 intercept 1
#'@param int3 intercept 2
#'@param beta1 slope 1
#'@param beta2 slope 2
#'@return a metabolic value

funtorp2 <- function(x, Tt, intr, intc, betat, betac, Ym) {
  out <- (ifelse(x<Tt,intr +betat*x,intc*exp(betac*x)))*Ym
  return(out)
}
#' funnorm
#'
#'fit the model for normotermic bats
#'@name funnorm
#'@param x a temperature
#'@param int1 intercept 1
#'@param beta1 slope 1
#'@return a metabolic value
funnorm2 <- function(x, inte, betat, Ym) {
  (inte +betat*x)*Ym
  return(out)
}

