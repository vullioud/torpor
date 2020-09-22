#' get predictions
#'
#' [tor_predict()] provides the predicted metabolic rate (MR) and 95CI bounds at
#' a given ambient temperature (Ta), in euthermic and/or torpid state.
#'
#'@name tor_predict
#'@aliases tor_predict
#'@family predict
#'@param tor_obj a fitted model from [tor_fit()]
#'@param Ta a vector of temperature
#'@return a data.frame
#'@export
#'@examples
#'\dontrun{
#'data(test_data2)
#'test_mod <- tor_fit(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'MTNZ = 1.49,
#'TLC = 28.88,
#'fitting_options = list(nc = 2, nb = 3000, ni = 5000))
#'tor_predict(tor_obj, Ta = 10:35)
#'}
tor_predict <- function(tor_obj, Ta){   #### En dessus de TLC = MTNZ, warning en dehors de MTNZ si plus que TLC.


  # retrieved the posterior of interest
  mod <- tor_obj$mod_parameter

  betat <- mod$sims.list$betat
  betac <- mod$sims.list$betac
  inte <- mod$sims.list$inte
  intc <- mod$sims.list$intc
  intr <- mod$sims.list$intr
  Tt <- mod$sims.list$Tt

  ##
  tlc <- tor_obj$out_mtnz_tlc$tlc_estimated
  Ym <- tor_obj$data$Ym


  X <- length(Ta)
  Ymean_t <- rep(NA, X)
  Y975_t <- rep(NA, X)
  Y025_t <- rep(NA, X)

  Ymean_n <- rep(NA, X)
  Y975_n <- rep(NA, X)
  Y025_n <- rep(NA, X)

  Ymean_b <- rep(NA, X)
  Y975_b <- rep(NA, X)
  Y025_b <- rep(NA, X)

  for (i in 1:X) {

    if(Ta[i] < tor_obj$out_mtnz_tlc$tlc_estimated) {


    Ymean_t[i] <- stats::median(tor_predict_fun(Ta[i], Tt, intr, intc, betat, betac, Ym))
    Y975_t[i] <- stats::quantile(tor_predict_fun(Ta[i],Tt, intr, intc, betat, betac, Ym),0.975)
    Y025_t[i] <- stats::quantile(tor_predict_fun(Ta[i],Tt, intr, intc, betat, betac, Ym),0.025)

    Ymean_n[i] <- stats::median(eut_predict_fun(Ta[i], inte, betat, Ym))
    Y975_n[i] <- stats::quantile(eut_predict_fun(Ta[i], inte, betat, Ym),0.975)
    Y025_n[i] <- stats::quantile(eut_predict_fun(Ta[i], inte, betat, Ym),0.025)
    } else {

      MTNZ <- tor_obj$out_mtnz_tlc$mtnz_points

      Ymean_b[i] <-  mean <- mean(MTNZ, na.rm = TRUE)
      Y975_b[i] <- mean(MTNZ)+1.96*stats::sd(MTNZ)/sqrt(length(tor_obj$out_mtnz_tlc$mtnz_points))
      Y025_b[i] <- mean(MTNZ)-1.96*stats::sd(MTNZ)/sqrt(length(tor_obj$out_mtnz_tlc$mtnz_points))

    }
  }


  #### values for torpor
  out_tor <- data.frame(Ta = Ta,
                     group = rep("Torpor", length(X)),
                     pred =  Ymean_t,
                     upr_95 = Y975_t,
                     lwr_95 = Y025_t)

  #### values for euthermy
  out_eut <- data.frame(Ta = Ta,
                     group = rep("Euthermy", length(X)),
                     pred =  Ymean_n,
                     upr_95 = Y975_n,
                     lwr_95 = Y025_n)


  ### values for > Tlc.
  out_mtnz <- data.frame(Ta = Ta,
                        group = rep("MTNZ", length(X)),
                        pred =  Ymean_b,
                        upr_95 = Y975_b,
                        lwr_95 = Y025_b)


  ###
  if(any(Ta > tor_obj$out_mtnz_tlc$tlc_estimated)) warning("Tuc is not considered: MTNZ is calculated independently of Ta above Tlc")

  if(!any(tor_obj$assignation$G == 1)){ ## no torpor
    out <- rbind(out_eut,out_mtnz)
  } else if (!any(tor_obj$assignation$G == 2)){ ## no euthermy
    out <- rbind(out_tor, out_mtnz)
  } else { ## both torpor and euthermy
    out <- rbind(out_tor, out_eut, out_mtnz)
  }

  return(stats::na.omit(out))

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
#'@param tor_obj a fitted model from [tor_fit()]
#'@return a data.frame
#'@export

tor_classify <- function(tor_obj){

data <- data.frame(measured_MR = tor_obj$data$Y,
                measured_Ta = tor_obj$data$Ta,
                predicted_state = tor_obj$assignation$G)

## sort tout MR et Ta et classification des data sur data step_4. 1 torpor,
## 2 euth, 3 NTMZ, 0 not classify. + Probastatus. + valeur predict qui vient de step_4.
mod <- tor_obj$mod_parameter

betat <- mod$sims.list$betat
betac <- mod$sims.list$betac
inte <- mod$sims.list$inte
intc <- mod$sims.list$intc
intr <- mod$sims.list$intr
Tt <- mod$sims.list$Tt
tlc <- mod$sims.list$tlc
Ym <- tor_obj$data$Ym

X <- nrow(data)
data$predicted_MR <- rep(NA, X)

data$classification <- dplyr::case_when(data$predicted_state == 0 ~ "Undefined",
                                     data$predicted_state == 1 ~ "Torpor",
                                     data$predicted_state == 2 ~ "Euthermy",
                                     data$predicted_state == 3 ~ "MTNZ")

for(i in 1:nrow(data)) {
  if (data$predicted_state[i] == 1) {
  data$predicted_MR[i] <- stats::median(tor_predict_fun(data$measured_Ta[i], Tt, intr, intc, betat, betac, Ym))
  } else if(data$predicted_state[i] == 2) {
  data$predicted_MR[i] <- stats::median(eut_predict_fun(data$measured_Ta[i], inte, betat, Ym))
  } else if(data$predicted_state[i] == 3) {
    data$predicted_MR[i] <- tor_obj$out_mtnz_tlc$mtnz_estimated
  } else {
    data$predicted_MR[i] <- NA
  }
}

data <- data[,c(1,2,5,4)]
return(data)
}
################################################################################


