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
#'@param tor_obj a fitted model from [tor_fit()]
#'@return a list of two data.frame. The first element returns the parameter estimates, the second reports the overlap values.
#'@export
#'@examples
#'data(test_data)
#'test2 <- tor_fit(Y = test_data[,2],
#'Ta = test_data[, 1],
#'fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'tor_summarise(test2)
tor_summarise <- function(tor_obj){

list(params = get_parameters(tor_obj),
     ppo = tor_ppo(tor_obj))## round print to 3 digit.
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
#'@param tor_obj a fitted model from [tor_fit()]
#'@return a data.frame
#'@export

tor_ppo <- function(tor_obj){
  overlap <- NULL
  ### TLC
  mod_tlc <- tor_obj$out_mtnz_tlc$model_1

  nbsamples <- nrow(mod_tlc$samples[[1]])*length(mod_tlc$samples)

  MIN<- max(tor_obj$out_mtnz_tlc$Ta2)
  MAX<-max(tor_obj$data$Ta)
  PR<- stats::runif(nbsamples,MIN,MAX)
  tlc_chain <- tor_obj$out_mtnz_tlc$tlc_distribution
  overlapTlc <- as.numeric(round(overlapping::overlap(x = list(tlc_chain, PR))$OV, digits = 3))
  out_tlc <- data.frame(name = "Tlc", overlap = overlapTlc)

  ## MR
  MIN <- 0
  MAX <- tor_obj$out_mtnz_tlc$mtnz_estimated
  PR<- stats::runif(nbsamples,MIN,MAX)
  MR_chain <- tor_obj$mod_parameter$sims.list$MR
  overlapMRr <- as.numeric(round(overlapping::overlap(x = list(MR_chain, PR))$OV, digits = 3))
  out_MR <- data.frame(name = "MR", overlap = overlapMRr)

  ## Tbe
  MIN <- tor_obj$out_mtnz_tlc$tlc_estimated
  MAX <- 100
  PR <- stats::runif(nbsamples,MIN,MAX)
  Tbe_chain <- tor_obj$mod_parameter$sims.list$Tbe
  overlapTbe <- as.numeric(round(overlapping::overlap(x = list(Tbe_chain, PR))$OV, digits = 3))

  out_Tbe <- data.frame(name = "Tbe", overlap = overlapTbe)

  ## TMR
  MIN <- 0
  for(i in 1:nbsamples){

    PR[i] <- stats::runif(1,MIN,tor_obj$mod_parameter$sims.list$MR[i])}

  TMR_chain <- tor_obj$mod_parameter$sims.list$TMR
  overlapTMR <- as.numeric(round(overlapping::overlap(x = list(TMR_chain, PR))$OV, digits = 3))
  out_TMR <- data.frame(name = "TMR", overlap = overlapTMR)

  ## tbt
  MIN<- -10
  MAX<- tor_obj$out_mtnz_tlc$tlc_estimated
  PR<- stats::runif(nbsamples,MIN,MAX)
  Tbt_chain <- tor_obj$mod_parameter$sims.list$Tbt
  overlapTbt <- as.numeric(round(overlapping::overlap(x = list(Tbt_chain, PR))$OV, digits = 3))
  out_Tbt <- data.frame(name = "Tbt", overlap = overlapTbt)

  out <- dplyr::bind_rows(out_tlc, out_MR, out_Tbe, out_TMR, out_Tbt) %>%
    dplyr::mutate(overlap = overlap *100) %>%
    dplyr::rename(ppo = overlap)
 ## add loop on out to flag overlap > 90 and 60 %.
  for (i in 1:nrow(out)){

    if (out$ppo[i] > 80) {
      warning(paste(out$name[i], "is not identifiable: 80% > PPO"))
    }
  }
  out
}
##############################################################################

#'Fetch parameters from model output
#'
#'[get_parameters()] Retrieves the parameters estimates from the model output
#'and creates a dataframe used in the [tor_summarise()] function
#'
#'@aliases get_parameters
#'@family summary
#'@param tor_obj a fitted model from [tor_fit()]
#'@return a data.frame
#'@export
get_parameters <- function(tor_obj){  ## out4 et out2 pour tlc.
  .data <- NULL
  Ym <- tor_obj$data$Ym

  ### first param for step_4
  mod_params <- tor_obj$mod_parameter

  params_mod_parameter <- c("tau", "inte", "intc","intr", "betat", "betac",
                            "Tt", "TMR", "MR", "tauy1", "tauy2", "Tbe", "Tbt") ## params of interest

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

  ##### add TLC and MTNZ
  mod_tlc <- tor_obj$out_mtnz_tlc$model_1

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

  ##### add MTNZ
  MTNZ <- tor_obj$out_mtnz_tlc$mtnz_points
  mean <- mean(MTNZ)
  CI_97.5 <- mean(MTNZ)+1.96*stats::sd(MTNZ)/sqrt(length(tor_obj$out_mtnz_tlc$mtnz_points))
  median <- stats::median(MTNZ)
  CI_2.5 <- mean(MTNZ)-1.96*stats::sd(MTNZ)/sqrt(length(tor_obj$out_mtnz_tlc$mtnz_points))
  Rhat <- NA

  x_mtnz <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
  x_mtnz$parameter <- "MTNZ" ## unify names
  x_mtnz <- x_mtnz[,c(6,1,2,3,4, 5)] ## put name as first column
  rownames(x_mtnz) <- NULL

  out <- rbind(x, x_tlc, x_mtnz)

  params_to_multiply <- c("tau1", "tau2", "inte", "intc", "intr", "MTNZ", "MR", "TMR", "betat")


  out <- out %>%
    dplyr::mutate(mean = ifelse(.data$parameter %in% params_to_multiply, .data$mean*Ym, .data$mean),
                  CI_2.5 = ifelse(.data$parameter %in% params_to_multiply, .data$CI_2.5*Ym, .data$CI_2.5),
                  median = ifelse(.data$parameter %in% params_to_multiply, .data$median*Ym, .data$median),
                  CI_97.5 = ifelse(.data$parameter %in% params_to_multiply, .data$CI_97.5*Ym, .data$CI_97.5)) %>%
    dplyr::mutate_if(is.numeric, ~ round(.x, digits = 3))

  warning("Values for MTNZ are directly computed from data points")
  return(out)
}


