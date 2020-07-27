#' Model summary
#'
#'[tor_summarise()] provides a comprehensive summary of the output
#'returned by [tor_fit()]. Values of mean, 95CI bounds, median and
#'Brooks–Gelman–Rubin criterion (i.e. chain convergence check) of parameters
#'posterior distributions are provided. Additionally, prior posterior overlap
#'values for parameters tlc, Tt and Betat are generated.
#'
#'@name tor_summarise
#'@aliases tor_summarise
#'@family summary
#'@param mod a fitted model from [tor_fit()]
#'@return a list of two data.frame. The first element returns the parameter estimates, the second reports the overlap values.
#'@export
#'@examples
#'data(test_data)
#'test2 <- estimate_parameters(Y = test_data[,2],Ta = test_data[, 1], fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'tor_summarise(tor_obj)
tor_summarise <- function(mod){

#out <- list()
## parameters
get_parameters(mod)  ## round print to 3 digit.
return(out)
}
##############################################################################

#' Check priors/posteriors overlap
#'
#'[tor_overlap()] generates prior/posterior overlap values for Tlc,
#'Tt and Betat. Values larger than 0.3 should lead to the conclusion that conforming
#'torpor, regulated torpor or thermoregulation respectively could not be modeled
#'with the data provided (Fasel et al. (in prep)). [tor_overlap()] can be used independently but is also used internally in
#'[tor_summarise()].
#'
#'@name tor_overlap
#'@aliases tor_overlap
#'@family summary
#'@param mod a fitted model from [tor_fit()]
#'@return a data.frame
#'@export
#'@examples
#'data(test_data2)
#'test <- tor_fit(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 1.49,
#'TLC = 28.8,
#'model = NULL,
#'fitting_options = list(nc = 1, ni = 8000, nb = 3000))
#'tor_overlap(mod = test)

tor_overlap <- function(mod){
# size <- nrow(mod$samples[[1]])
#
# tlc <- mod$mean$tlc
# bmr <- mod$data$BMR
# Ym  <- mod$data$Ym
#
# PRTintc <- truncnorm::rtruncnorm(size,a = 0, b = bmr/Ym, mean = 0, sd = sqrt(1000))  ## variable based on model flexible
# PRTt <- truncnorm::rtruncnorm(size, b = tlc, mean = 0, sd = sqrt(1000))
# PRbeta <- truncnorm::rtruncnorm(size, b = 0, mean = 0, sd = sqrt(100))
#
# intc_overlap <- get_overlap(mod, params = "intc", priors = PRTintc)
# Tt_overlap <- get_overlap(mod, params = "Tt", priors = PRTt)
# betat_overlap <- get_overlap(mod, params = "betat", priors = PRbeta)

#
#   prior_tlc <- truncnorm::rtruncnorm(20000,a = mod$mean$TLC,mean=0,sd=100)
#   prior_Tt <- truncnorm::rtruncnorm(20000,b = mod$mean$TLC,mean=0,sd=100)
#   prior_betat <- truncnorm::rtruncnorm(20000,b = 0,mean=0,sd=100)
#
#   tlc_overlap <- get_overlap(mod, "tlc", prior_tlc)
#   Tt_overlap <- get_overlap(mod, "Tt", prior_Tt)
#   betat_overlap <- get_overlap(mod, "betat", prior_betat)

  # if (intc_overlap >= 0.3) {
  #   warning("Parameters intc, betac, and intc not identifiable")
  # }
  # if (Tt_overlap >= 0.3) {
  #   warning("Parameters Tt and intr not identifiable")
  # }
  # if (betat_overlap >= 0.3) {
  #   warning("Parameters tlc, betac, and intc not identifiable")
  # }
  #
  # out <- data.frame(parameter = c("intc", "Tt", "Betat"),
  #                   overlap = c(intc_overlap, Tt_overlap, betat_overlap))
  #
  # return(out)
}
##############################################################################

#'Fetch parameters from model output
#'
#'[get_parameters()] Retrieves the parameters estimates from the model output
#'and creates a dataframe used in the [tor_summarise()] function
#'
#'@aliases get_parameters
#'@family summary
#'@param mod a fitted model from [tor_fit()]
#'@return a data.frame
#'@export
#'@examples
#'data(test_data2)
#'test <- tor_fit(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 1.49,
#'TLC = 28.8,
#'model = NULL,
#'fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'get_parameters(mod = test)
get_parameters <- function(tor_obj) {  ## out4 et out2 pour tlc.

  Ym <- tor_obj$data$Ym

  ### first param for step_4
  mod_params <- tor_obj$mod_parameter

  params_mod_parameter <- c("tau", "inte", "intc", "intr", "betat", "betac", "Tt", "TMR", "MRr", "coef") ## params of interest

  mean <- unlist(mod_params$mean[params_mod_parameter])
  CI_97.5 <- unlist(mod_params$q97.5[params_mod_parameter])
  median <- unlist(mod_params$q50[params_mod_parameter])
  CI_2.5 <- unlist(mod_params$q2.5[params_mod_parameter])
  Rhat <- unlist(mod_params$Rhat[params_mod_parameter])

  ## frame the output in a df
  x <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
  x$parameter <- rownames(x)
  x <- x[,c(6,1,2,3,4, 5)] ## put name as first column
  rownames(x) <- NULL
  ##### add TLC and BMR
  mod_tlc <- tor_obj$out_bmr_tlc$model_1

  mean <- unlist(mod_tlc$mean["tlc"])
  CI_97.5 <- unlist(mod_tlc$q97.5["tlc"])
  median <- unlist(mod_tlc$q50["tlc"])
  CI_2.5 <- unlist(mod_tlc$q2.5["tlc"])
  Rhat <- unlist(mod_tlc$Rhat["tlc"])

  x_tlc <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
  x_tlc$parameter <- "Tlc" ## unify names
  x_tlc <- x_tlc[,c(6,1,2,3,4, 5)] ## put name as first column
  rownames(x_tlc) <- NULL

  out <- rbind(x, x_tlc)

  ##### add bmr
  bmr <- tor_obj$out_bmr_tlc$bmr_points
  mean <- mean(bmr)
  CI_97.5 <- mean(bmr)+1.96*sd(bmr)/sqrt(length(tor_obj$out_bmr_tlc$bmr_points))
  median <- median(bmr)
  CI_2.5 <- mean(bmr)-1.96*sd(bmr)/sqrt(length(tor_obj$out_bmr_tlc$bmr_points))
  Rhat <- NA

  x_bmr <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
  x_bmr$parameter <- "Bmr" ## unify names
  x_bmr <- x_bmr[,c(6,1,2,3,4, 5)] ## put name as first column
  rownames(x_bmr) <- NULL

  out <- rbind(x, x_tlc, x_bmr)

  params_to_multiply <- c("tau1", "tau2", "inte", "intc", "intr", "BMR", "MRr", "TMR", "betat")

  out <- out %>%
    dplyr::mutate(mean = ifelse(.data$parameter %in% params_to_multiply, .data$mean*Ym, .data$mean),
                  CI_2.5 = ifelse(.data$parameter %in% params_to_multiply, .data$CI_2.5*Ym, .data$CI_2.5),
                  median = ifelse(.data$parameter %in% params_to_multiply, .data$median*Ym, .data$median),
                  CI_97.5 = ifelse(.data$parameter %in% params_to_multiply, .data$CI_97.5*Ym, .data$CI_97.5))


  warning("Values for MTNZ are directly computed from data points")
  return(out)
}
##############################################################################

#' get_overlap // internal
#'
#' Compare the distribution overlap between a prior and a posterior // used internally only
#'
#'@name get_overlap
#'@aliases get_overlap
#'@family summary
#'@param mod a fitted model from [tor_fit()]
#'@param params a parameter of the model
#'@param priors a posterior distribution (as. numerical vector)
#'@importFrom magrittr %>%
#'@return a numeric value
get_overlap <- function(mod, params, priors) {
#
#   chains <- purrr::map_dfr(mod$samples, ~ as.data.frame(.x) %>%
#                              dplyr::select(paste(params)))
#
#   as.numeric(round(overlapping::overlap(x = list(chains[,1], priors))$OV, digits = 3))
}
##############################################################################
