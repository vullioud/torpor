#' Model predictions
#'
#' [tor_predict()] The function provides the predicted M and 95% credible interval
#' boundaries at a defined Ta given a certain model, in euthermic and/or torpid state.
#'
#' Note that the predictions are computed slightly differently depending on whether one uses `CI = TRUE` or `CI = FALSE`.
#' In the former case (the default), for each value of Ta, one prediction is computed for every single iteration of the posteriors and the median value of these predictions is then computed.
#' In the latter case, for each value of Ta, a single prediction is directly computed by assuming that parameters equate their median values from the corresponding posterior distribution.
#'
#'@name tor_predict
#'@aliases tor_predict
#'@family predict
#'@param tor_obj a fitted model from [tor_fit()]
#'@param Ta a vector of temperature
#'@param CI whether or not to compute prediction intervals (default = TRUE, which reduces the speed of computation noticeably for large jobs)
#'@return a data.frame
#'@export
#'@examples
#'\dontrun{
#' test_mod <- tor_fit(M = test_data2[,2],
#'                     Ta = test_data2[, 1],
#'                     fitting_options = list(nc = 1, nb = 3000, ni = 5000))
#' tor_predict(test_mod, Ta = 10:35)
#' tor_predict(test_mod, Ta = 10:35, CI = FALSE)
#'}
tor_predict <- function(tor_obj, Ta, CI = TRUE){
  if (!("tor_obj" %in% class(tor_obj))) stop("tor_obj need to be of class tor_obj")
  if (!(class(Ta) %in% c("integer", "numeric"))) stop("Ta need to be a integer or a numeric vector")

  # retrieved the posterior of interest
  mod <- tor_obj$mod_parameter

  betat <- mod$sims.list$betat
  betac <- mod$sims.list$betac
  inte <- mod$sims.list$inte
  intc <- mod$sims.list$intc
  intr <- mod$sims.list$intr
  Tt <- mod$sims.list$Tt

  ##
  Tlc <- tor_obj$out_Mtnz_Tlc$Tlc_estimated
  Ym <- tor_obj$data$Ym
  Mtnz <- tor_obj$out_Mtnz_Tlc$Mtnz_estimated

  X <- length(Ta)

  if (!CI) {
    # vectorial computation for simple case:
    Ymean_t <- tor_predict_fun(Ta, stats::median(Tt), stats::median(intr), stats::median(intc), stats::median(betat), stats::median(betac), Ym)
    Ymean_n <- eut_predict_fun(Ta, stats::median(inte), stats::median(betat), Ym)
    Ymean_b <- mean(Mtnz, na.rm = TRUE)
    Y975_t <- Y025_t <- Y975_n <- Y025_n <- Y975_b <- Y025_b <- NA

    } else {
    # non-vectorial computation:
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

      if (Ta[i] < Tlc) {


      Ymean_t[i] <- stats::median(tor_predict_fun(Ta[i], Tt, intr, intc, betat, betac, Ym))
      Y975_t[i] <- stats::quantile(tor_predict_fun(Ta[i],Tt, intr, intc, betat, betac, Ym), 0.975)
      Y025_t[i] <- stats::quantile(tor_predict_fun(Ta[i],Tt, intr, intc, betat, betac, Ym), 0.025)

      Ymean_n[i] <- stats::median(eut_predict_fun(Ta[i], inte, betat, Ym))
      Y975_n[i] <- stats::quantile(eut_predict_fun(Ta[i], inte, betat, Ym),0.975)
      Y025_n[i] <- stats::quantile(eut_predict_fun(Ta[i], inte, betat, Ym),0.025)

      } else {

        if (length(Mtnz < 2)) { # cases when Mtnz is provided by the user
          Ymean_b[i] <- Mtnz
          Y975_b[i] <- Mtnz
          Y025_b[i] <- Mtnz

        } else {
        Ymean_b[i] <- mean(Mtnz, na.rm = TRUE)
        Y975_b[i] <- mean(Mtnz) + 1.96*stats::sd(Mtnz)/sqrt(length(tor_obj$out_Mtnz_Tlc$Mtnz_points))
        Y025_b[i] <- mean(Mtnz) - 1.96*stats::sd(Mtnz)/sqrt(length(tor_obj$out_Mtnz_Tlc$Mtnz_points))
    }
      }
    }
  }

  #### values for torpor
  out_tor <- data.frame(Ta = Ta,
                     assignment = rep("Torpor", length(X)),
                     pred =  Ymean_t,
                     upr_95 = Y975_t,
                     lwr_95 = Y025_t)

  #### values for euthermia
  out_eut <- data.frame(Ta = Ta,
                     assignment = rep("Euthermia", length(X)),
                     pred =  Ymean_n,
                     upr_95 = Y975_n,
                     lwr_95 = Y025_n)


  ### values for > Tlc.
  out_Mtnz <- data.frame(Ta = Ta,
                        assignment = rep("Mtnz", length(X)),
                        pred =  Ymean_b,
                        upr_95 = Y975_b,
                        lwr_95 = Y025_b)
  if (!CI) {
    # specific filtering needed since we create Mtnz unconditionally in the simple case
    out_Mtnz <- out_Mtnz[out_Mtnz$Ta >= Tlc, ]
  }


  ###
  if (any(Ta > tor_obj$out_Mtnz_Tlc$Tlc_estimated)) warning("Tuc is not considered: Mtnz is calculated independently of Ta above Tlc")

  if (!any(tor_obj$assignation$G == 1)) { ## no torpor
    out <- rbind(out_eut, out_Mtnz)
  } else if (!any(tor_obj$assignation$G == 2)) { ## no euthermia
    out <- rbind(out_tor, out_Mtnz)
  } else { ## both torpor and euthermia
    out <- rbind(out_tor, out_eut, out_Mtnz)
  }

  out[!is.na(out$pred), ]

}
################################################################################
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
#'function to fit the model in euthermia
#'@name eut_predict_fun
#'@family predict
#'@param x a temperature
#'@param inte intercept 1
#'@param betat slope 1
#'@param Ym mean of Y to back transform
#'@return a numerical value
eut_predict_fun <- function(x, inte, betat, Ym) { ## backtransform parameters to prediction for euthermia
  out <- (inte +betat*x)*Ym
  return(out)
}
################################################################################

#' Assign metabolic state
#'
#'[tor_assign()] The function assign the individual points according to their
#'estimated state.
#'
#'@name tor_assign
#'@aliases tor_assign
#'@family predict
#'@param tor_obj a fitted model from [tor_fit()]
#'@return a data.frame
#'@export

tor_assign <- function(tor_obj){
  if(!("tor_obj" %in% class(tor_obj))) stop("tor_obj need to be of class tor_obj")

data <- data.frame(measured_M = tor_obj$data$Y,
                measured_Ta = tor_obj$data$Ta,
                predicted_state = tor_obj$assignation$G)

## sort tout M et Ta et classification des data sur data step_4. 1 torpor,
## 2 euth, 3 NTMZ, 0 not classify. + Probastatus. + valeur predict qui vient de step_4.
mod <- tor_obj$mod_parameter

betat <- mod$sims.list$betat
betac <- mod$sims.list$betac
inte <- mod$sims.list$inte
intc <- mod$sims.list$intc
intr <- mod$sims.list$intr
Tt <- mod$sims.list$Tt
Tlc <- mod$sims.list$Tlc
Ym <- tor_obj$data$Ym

X <- nrow(data)
data$predicted_M <- rep(NA, X)

data$assignment <- dplyr::case_when(data$predicted_state == 0 ~ "Undefined",
                                     data$predicted_state == 1 ~ "Torpor",
                                     data$predicted_state == 2 ~ "Euthermia",
                                     data$predicted_state == 3 ~ "Mtnz")

for(i in 1:nrow(data)) {
  if (data$predicted_state[i] == 1) {
  data$predicted_M[i] <- stats::median(tor_predict_fun(data$measured_Ta[i], Tt, intr, intc, betat, betac, Ym))
  } else if(data$predicted_state[i] == 2) {
  data$predicted_M[i] <- stats::median(eut_predict_fun(data$measured_Ta[i], inte, betat, Ym))
  } else if(data$predicted_state[i] == 3) {
    data$predicted_M[i] <- tor_obj$out_Mtnz_Tlc$Mtnz_estimated
  } else {
    data$predicted_M[i] <- NA
  }
}

data[,c(1,2,5,4)]
}


