#' tor_summarise
#'
#'[tor_summarise()] provides a comprehensive summary of the output
#'returned by [tor_fit()]. Values of mean, 95CI bounds, median and
#'Brooks–Gelman–Rubin criterion (i.e. chain convergence check) of parameters
#'posterior distributions are provided. Additionally, prior posterior overlap
#'values for different parameters are generated. Finally, the mean MR value used
#'to standardize the MR variable is reported.
#'
#'@name tor_summarise
#'@aliases tor_summarise
#'@family summary
#'@param mod a fitted model with [tor_fit()]
#'@return a list of 2. The first element is a tibble with the parameter estimates, the second is the overlap values
#'@export
#'@examples
#'data(test_data)
#'test2 <- tor_fit(MR = test_data[,2],
#'Ta = test_data[, 1],
#'BMR = 29,
#'TLC = 28.8,
#'model = NULL,
#'fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'tor_summarise(mod = test2)

tor_summarise <- function(mod){

out <- list()
## parameters
out$parameter_estimates <- get_parameters(mod)
## overlap
out$overlap <- tor_overlap(mod)
return(out)
}
##############################################################################

#' tor_overlap
#'
#'[tor_overlap()] generates prior/posterior overlap values for Tm,
#'Tt and Betat.
#'Note: Values larger than 0.3 should lead to the conclusion that conforming
#'torpor, regulated torpor or thermoregulation respectively could not be modeled
#' with the data provided (Fasel et al. (in prep)). Can be used independently but is also used internally in
#' tor_summarise()
#'@name tor_overlap
#'@aliases tor_overlap
#'@family summary
#'@param mod a fitted model with [tor_fit()]
#'@return a list of 3. With overlapping values for tlc, Tt and Betat
#'@export
#'@examples

#'data(test_data2)
#'test <- tor_fit(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 1.49,
#'TLC = 28.8,
#'model = NULL,
#'fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'tor_overlap(mod = test)

tor_overlap <- function(mod){

  prior_tlc <- truncnorm::rtruncnorm(20000,a= mod$mean$TLC,mean=0,sd=100)
  prior_Tt <- truncnorm::rtruncnorm(20000,b=mod$mean$TLC,mean=0,sd=100)
  prior_betat <- truncnorm::rtruncnorm(20000,b=0,mean=0,sd=100)

  tlc_overlap <- get_overlap(mod, "tlc", prior_tlc)
  Tt_overlap <- get_overlap(mod, "Tt", prior_Tt)
  betat_overlap <- get_overlap(mod, "betat", prior_betat)

  if(tlc_overlap >= 0.3){
    warning("Parameters tlc, betac, and intc not identifiable")
  }
  if(Tt_overlap >= 0.3){
    warning("Parameters Tt and intr not identifiable")
  }
  if(betat_overlap >= 0.3){
    warning("Parameters tlc, betac, and intc not identifiable")
  }

  out <- data.frame(parameter = c("tlc", "Tt", "Betat"),
                    overlap = c(tlc_overlap, Tt_overlap, betat_overlap))

  return(out)
}
##############################################################################

#'get_parameters // internal
#'
#'Retrieves the parameters estimates from the model output
#'it creates a dataframe used in the [tor_summarise()] function
#'@aliases get_parameters
#'@family summary
#'@param mod a fitted model with [tor_fit()]
#'@return a data.frame


get_parameters <- function(mod) {

  params <- c("tauy", "inte", "intc", "intr", "betat", "betac", "Tt", "tlc") ## params of interest

  mean <- unlist(mod$mean[params])
  CI_97.5 <- unlist(mod$q97.5[params])
  median <- unlist(mod$q50[params])
  CI_2.5 <- unlist(mod$q2.5[params])
  Rhat <- unlist(mod$Rhat[params])

  ## frame the output in a df
  x <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
  x$parameter <- rownames(x)
  x$parameter[x$parameter == "tlc"] <- "Tm" ## unify names
  x <- x[,c(6,1,2,3,4, 5)] ## put name as first column
  rownames(x) <- NULL
  return(x)
}
##############################################################################

#' get_overlap // internal
#'
#' Compare the distribution overlap between a prior and a posterior // used internally only
#'@name get_overlap
#'@aliases get_overlap
#'@family summary
#'@param mod a fitted model with [tor_fit()]
#'@param params a parameter of the model
#'@param priors a posterior distribution (as. numerical vector)
#'@importFrom magrittr %>%
#'@return a numeric value
get_overlap <- function(mod, params, priors) {

  chains <- purrr::map_dfr(mod$samples, ~ as.data.frame(.x) %>%
                             dplyr::select(paste(params)))

  as.numeric(round(overlapping::overlap(x = list(chains[,1], priors))$OV, digits = 3))
}
##############################################################################
