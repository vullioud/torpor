#' get predictions
#'
#' [tor_predict()] provides the predicted metabolic rate (MR) and 95CI bounds at
#' a given ambient temperature (Ta), in euthermic and/or torpid state.
#'
#'@name tor_predict
#'@aliases tor_predict
#'@family predict
#'@param mod a fitted model from [tor_fit()]
#'@param Ta a vector of temperature
#'@return a data.frame
#'@export
#'@examples
#'\dontrun{
#'data(test_data2)
#'test_mod <- tor_fit(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 1.49,
#'TLC = 28.88,
#'fitting_options = list(nc = 2, nb = 3000, ni = 5000))
#'tor_predict(mod = test_mod, Ta = 1:20)
#'}
tor_predict <- function(mod, Ta){


  # retrieved the posterior of interest
  betat <- mod$sims.list$betat
  betac <- mod$sims.list$betac
  inte <- mod$sims.list$inte
  intc <- mod$sims.list$intc
  intr <- mod$sims.list$intr
  Tt <- mod$sims.list$Tt
  tlc <- mod$sims.list$tlc
  Ym <- mod$data$Ym


  X <- length(Ta)
  Ymeant <- rep(NA, X)
  Y975t <- rep(NA, X)
  Y025t <- rep(NA, X)

  Ymeann <- rep(NA, X)
  Y975n <- rep(NA, X)
  Y025n <- rep(NA, X)

  for (i in 1:X) {
    Ymeant[i] <- stats::median(tor_predict_fun(Ta[i], Tt, intr, intc, betat, betac, Ym))
    Y975t[i] <- stats::quantile(tor_predict_fun(Ta[i],Tt, intr, intc, betat, betac, Ym),0.975)
    Y025t[i] <- stats::quantile(tor_predict_fun(Ta[i],Tt, intr, intc, betat, betac, Ym),0.025)
    Ymeann[i] <- stats::median(eut_predict_fun(Ta[i], inte, betat, Ym))
    Y975n[i] <- stats::quantile(eut_predict_fun(Ta[i], inte, betat, Ym),0.975)
    Y025n[i] <- stats::quantile(eut_predict_fun(Ta[i], inte, betat, Ym),0.025)
  }


  #### values for torpor
  out_tor <- data.frame(Ta = Ta,
                     group = rep("Torpor", length(X)),
                     pred =  Ymeant,
                     upr_95 = Y975t,
                     lwr_95 = Y025t)

  #### values for euthermy
  out_eut <- data.frame(Ta = Ta,
                     group = rep("Euthermy", length(X)),
                     pred =  Ymeann,
                     upr_95 = Y975n,
                     lwr_95 = Y025n)



  if(length(mod$mean$G[mod$mean$G > 1.5]) == 0){ ## no torpor
    out <- out_eut
  } else if (length(mod$mean$G[mod$mean$G < 1.5]) == 0){ ## no euthermy
    out <- out_tor
  } else { ## both torpor and euthermy
    out <- rbind(out_tor, out_eut)
  }

  return(out)

}
################################################################################
#'
#' tor_predict_fun // internal
#'
#'Function to fit the model in torpor
#'@name tor_predict_fun
#'@family predict
#'@param x a temperature
#'@param Tt turning point T
#'@param intr intercept 1
#'@param intc intercept 2
#'@param betat slope 1
#'@param betac slope 2
#'@param Ym mean of Y to back-transform
#'@return a numerical value

tor_predict_fun <- function(x, Tt, intr, intc, betat, betac, Ym) {
  out <- (ifelse(x<Tt,intr +betat*x,intc*exp(betac*x)))*Ym ## backtransform parameters to prediction for torpor
  return(out)
}
################################################################################


#' eut_predict_fun // internal
#'
#'function to fit the model in euthermy
#'@name eut_predict_fun
#'@family predict
#'@param x a temperature
#'@param inte intercept 1
#'@param betat slope 1
#'@param Ym mean of Y to back transform
#'@return a numerical value
eut_predict_fun <- function(x, inte, betat, Ym) { ## backtransform parameters to prediction for euthermy
  out <- (inte +betat*x)*Ym
  return(out)
}
################################################################################


#' tor_classify
#'
#'[tor_classify()] returns the raw data with the related
#'predicted state values (between 1 and 2), which leads to the state
#'classification (torpor or euthermy). Additionally, it also provides the
#'predicted metabolic rate (MR) at the given ambient temperature (Ta).
#'
#'@name tor_classify
#'@aliases tor_classify
#'@family predict
#'@param mod a fitted model from [tor_fit()]
#'@return a data.frame
#'@export
#'@examples
#'data(test_data2)
#'test_mod <- tor_fit(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 1.49,
#'TLC = 28.88,
#'fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'tor_classify(mod = test_mod)

tor_classify <- function(mod){
x <- data.frame(measured_MR = (mod$data$Y)*mod$data$Ym[1],
                measured_Ta = mod$data$Ta,
                predicted_state = mod$mean$G)

betat <- mod$sims.list$betat
betac <- mod$sims.list$betac
inte <- mod$sims.list$inte
intc <- mod$sims.list$intc
intr <- mod$sims.list$intr
Tt <- mod$sims.list$Tt
tlc <- mod$sims.list$tlc
Ym <- mod$data$Ym


X <- nrow(x)
x$predicted_MR <- rep(NA, X)
x$classification <- ifelse(x$predicted_state > 1.5, "Torpor", "Euthermy")

for(i in 1:nrow(x)) {
  if (x$predicted_state[i] > 1.5) {
  x$predicted_MR[i] <- stats::median(tor_predict_fun(x$measured_Ta[i], Tt, intr, intc, betat, betac, Ym))
  } else {
  x$predicted_MR[i] <- stats::median(eut_predict_fun(x$measured_Ta[i], inte, betat, Ym))
  }
}
return(x)
}
################################################################################


