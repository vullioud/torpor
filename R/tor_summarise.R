#' Model summary
#'
#'[tor_summarise()] provides a comprehensive summary of the output
#'returned by [tor_fit()]. Values of mean, 95CI bounds, median and
#'Brooks–Gelman–Rubin criterion (i.e. chain convergence check) of parameters
#'posterior distributions are provided. Additionally, prior posterior overlap
#'values for parameters Tlc, Tt and Betat are generated.
#'
#'@name tor_summarise
#'@aliases tor_summarise
#'@family summary
#'@param tor_obj a fitted model from [tor_fit()]
#'@return a list of two data.frame. The first element returns the parameter estimates, the second reports the overlap values.
#'@export
#'@examples
#'test2 <- tor_fit(M = test_data[,2],
#'                 Ta = test_data[, 1],
#'                 fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'tor_summarise(test2)
tor_summarise <- function(tor_obj){

  if(!("tor_obj" %in% class(tor_obj))) stop("tor_obj need to be of class tor_obj")

  list(params = get_parameters(tor_obj),
       ppo = tor_ppo(tor_obj))## round print to 3 digit.
}
##############################################################################

#' Check priors/posteriors overlap
#'
#'[tor_ppo()] generates prior/posterior overlap values for Tlc,
#'Tt and Betat. Values larger than 0.3 should lead to the conclusion that conforming
#'torpor, regulated torpor or thermoregulation respectively could not be modeled
#'with the data provided (Fasel et al. (in prep)). [tor_ppo()] can be used independently but is also used internally in
#'[tor_summarise()].
#'
#'@name tor_ppo
#'@aliases tor_ppo
#'@family summary
#'@param tor_obj a fitted model from [tor_fit()]
#'@return a data.frame
#'@export

tor_ppo <- function(tor_obj){

  if(!("tor_obj" %in% class(tor_obj))) stop("tor_obj need to be of class tor_obj")

  ### Tlc
  mod_param <- tor_obj$mod_parameter
  mod_Tlc <- tor_obj$out_Mtnz_Tlc$model_1

  nbsamples <- nrow(mod_param$samples[[1]])*length(mod_param$samples)

  if(!(is.null(mod_Tlc))) {

    MIN<- max(tor_obj$out_Mtnz_Tlc$Ta2)
    MAX<-max(tor_obj$data$Ta)
    PR<- truncnorm::rtruncnorm(nbsamples,
                               a=MIN,
                               b=MAX,
                               mean=0,
                               sd=sqrt(1/0.001))
    Tlc_chain <- tor_obj$out_Mtnz_Tlc$Tlc_distribution
    overlapTlc <- as.numeric(round(overlapping::overlap(x = list(Tlc_chain, PR))$OV, digits = 3))
    out_Tlc <- data.frame(name = "Tlc", overlap = overlapTlc)

  } else {
    out_Tlc <- data.frame(name = "Tlc", overlap = 0)
  }

  ## MR
  PR<- rep(NA, nbsamples)
  MAX<-tor_obj$out_Mtnz_Tlc$Mtnz_estimated/tor_obj$data$Ym
  for(i in 1:nbsamples){
    PR[i] <- truncnorm::rtruncnorm(1,
                                   a=tor_obj$mod_parameter$sims.list$TMR[i],
                                   b=MAX,
                                   mean=0,
                                   sd=sqrt(1/0.001))
  }
  MR_chain <- tor_obj$mod_parameter$sims.list$MR
  overlapMR <- as.numeric(round(overlapping::overlap(x = list(MR_chain, PR))$OV, digits = 3))
  out_MR <- data.frame(name = "MR", overlap = overlapMR)

  ## Tbe
  MIN <- tor_obj$out_Mtnz_Tlc$Tlc_estimated
  MAX <- 50
  PR <- truncnorm::rtruncnorm(nbsamples,
                              a=MIN,
                              b=MAX,
                              mean=0,
                              sd=sqrt(1/0.001))
  Tbe_chain <- tor_obj$mod_parameter$sims.list$Tbe
  overlapTbe <- as.numeric(round(overlapping::overlap(x = list(Tbe_chain, PR))$OV, digits = 3))

  out_Tbe <- data.frame(name = "Tbe", overlap = overlapTbe)

  ## TMR
  MIN <- 0
  MAX <- 0.8 * tor_obj$out_Mtnz_Tlc$Mtnz_estimated/tor_obj$data$Ym
  PR<- truncnorm::rtruncnorm(nbsamples,
                             a=MIN,
                             b=MAX,
                             mean=0,
                             sd=sqrt(1/0.001))
  TMR_chain <- tor_obj$mod_parameter$sims.list$TMR
  overlapTMR <- as.numeric(round(overlapping::overlap(x = list(TMR_chain, PR))$OV, digits = 3))
  out_TMR <- data.frame(name = "TMR", overlap = overlapTMR)

  ## tbt
  MIN<- -5
  for(i in 1:nbsamples){
    PR[i]<- truncnorm::rtruncnorm(1,
                                  a=-5,
                                  mean=0,
                                  sd=sqrt(1/0.001))
  }

  Tbt_chain <- tor_obj$mod_parameter$sims.list$Tbt
  overlapTbt <- as.numeric(round(overlapping::overlap(x = list(Tbt_chain, PR))$OV, digits = 3))
  out_Tbt <- data.frame(name = "Tbt", overlap = overlapTbt)

  out <- dplyr::bind_rows(out_Tlc, out_MR, out_Tbe, out_TMR, out_Tbt) %>%
    dplyr::mutate(overlap = .data$overlap *100) %>%
    dplyr::rename(ppo = .data$overlap)
  ## add loop on out to flag overlap > 80 %.
  for (i in 1:nrow(out)){

    if (out$ppo[i] > 80) {
      warning(paste(out$name[i], "is not identifiable: PPO > 80%"))
    }
  }
  out
}
##############################################################################

#'Fetch parameters from model output
#'
#'[get_parameters()] Retrieves the parameters estimates from the model output
#'and creates a data.frame used in the [tor_summarise()] function
#'
#'@aliases get_parameters
#'@family summary
#'@param tor_obj a fitted model from [tor_fit()]
#'@return a data.frame
#'@importFrom rlang .data
#'@export
get_parameters <- function(tor_obj){  ## out4 et out2 pour Tlc.

  if(!("tor_obj" %in% class(tor_obj))) stop("tor_obj needs to be of class tor_obj")

  Ym <- tor_obj$data$Ym

  ### first param for step_4
  mod_params <- tor_obj$mod_parameter

  params_mod_parameter <- c("tau", "inte", "intc","intr", "betat", "betac",
                            "Tt", "TMR", "MR", "Tbe", "Tbt") ## params of interest

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

  ##### add Tlc and Mtnz
  mod_Tlc <- tor_obj$out_Mtnz_Tlc$model_1

  mean <- unlist(mod_Tlc$mean["Tlc"])
  CI_97.5 <- unlist(mod_Tlc$q97.5["Tlc"])
  median <- unlist(mod_Tlc$q50["Tlc"])
  CI_2.5 <- unlist(mod_Tlc$q2.5["Tlc"])
  Rhat <- unlist(mod_Tlc$Rhat["Tlc"])


  if(is.null(mod_Tlc)) {
    message("Tlc was not estimated from the data!")
    mean <- tor_obj$out_Mtnz_Tlc$Tlc_estimated
    CI_97.5 <- NA
    median <- tor_obj$out_Mtnz_Tlc$Tlc_estimated
    CI_2.5 <- NA
    Rhat <- NA
  }

  x_Tlc <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
  x_Tlc$parameter <- "Tlc" ## unify names
  x_Tlc <- x_Tlc[,c(6,1,2,3,4, 5)] ## put name as first column
  rownames(x_Tlc) <- NULL

  out <- rbind(x, x_Tlc)

  ##### add Mtnz
  Mtnz <- tor_obj$out_Mtnz_Tlc$Mtnz_points
  mean <- mean(Mtnz)
  CI_97.5 <- mean(Mtnz)+1.96*stats::sd(Mtnz)/sqrt(length(tor_obj$out_Mtnz_Tlc$Mtnz_points))
  median <- stats::median(Mtnz)
  CI_2.5 <- mean(Mtnz)-1.96*stats::sd(Mtnz)/sqrt(length(tor_obj$out_Mtnz_Tlc$Mtnz_points))
  Rhat <- NA

  x_Mtnz <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
  x_Mtnz$parameter <- "Mtnz" ## unify names
  x_Mtnz <- x_Mtnz[,c(6,1,2,3,4, 5)] ## put name as first column
  rownames(x_Mtnz) <- NULL

  out <- rbind(x, x_Tlc, x_Mtnz)

  params_to_multiply <- c("tau1", "tau2", "inte", "intc", "intr", "MR", "TMR", "betat")

  out <- out %>%
    dplyr::mutate(mean = ifelse(.data$parameter %in% params_to_multiply, .data$mean*Ym, .data$mean),
                  CI_2.5 = ifelse(.data$parameter %in% params_to_multiply, .data$CI_2.5*Ym, .data$CI_2.5),
                  median = ifelse(.data$parameter %in% params_to_multiply, .data$median*Ym, .data$median),
                  CI_97.5 = ifelse(.data$parameter %in% params_to_multiply, .data$CI_97.5*Ym, .data$CI_97.5)) %>%
    dplyr::mutate_if(is.numeric, ~ round(.x, digits = 3))

  out
}


