##############################################################################

#'complete_args // internal
#'
#'this function is used internally for tor_fit
#'From Alex Courtiol - informed.
#'
#'@name complete_args
#'@aliases complete_args
#'@param fn a function

.complete_args <- function(fn) {
  ## This function should not be called by the user.
  ## It keeps the default list elements when
  ## a new list with fewer elements is provided.
  env <- parent.frame()
  args <- formals(fn)
  for (arg_name in names(args)) {
    if (is.call(arg <- args[[arg_name]])) {
      if (arg[1] == "list()") {
        arg_input <- mget(names(args), envir = env)[[arg_name]]
        arg_full  <- eval(formals(fn)[[arg_name]])
        if (is.null(names(arg_input))) {
          if (length(arg_input) == length(arg_full)) {
            names(arg_input) <- names(arg_full)
          }
          else {
            stop(paste("The list", arg_name, "should contain names, or be of length equal to the default."))
          }
        }
        arg_full_updated <- utils::modifyList(arg_full, arg_input)
        assign(arg_name, arg_full_updated, envir = env)
      }
    }
  }
  return(NULL)
}


##############################################################################

#'Find the low Tlc.
#'
#' This function estimate the lowest Tlc and the Mtnz
#'
#'@name find_low_Tlc_Mtnz
#'@param M A vector of metabolic rate measure.
#'@param Ta A vector of ambient temperature.
#'@return A named vector with Mtmz and low_Tlc.
#'@export
#'@examples
#'t_1 <- find_low_Tlc_Mtnz(M = test_data$VO2, Ta = test_data$Ta)
find_low_Tlc_Mtnz <- function(M,Ta){
  ### clean arguments function
  Y <- M
  x <- stats::na.omit(cbind(Ta, Y))

  Y <- x[,2]
  Ta <- x[,1]

  SEQX <- seq(from = max(Ta, na.rm = TRUE), to = min(Ta, na.rm = TRUE), length.out = 100)

  Nb <- rep(NA, 100)
  whitetest <- rep(NA, 100)
  p <- rep(NA, 100)
  p2 <- rep(NA, 100)
  #Define low_Tlc and slope
  for (j in 1:100) {

    Nb[j] <- length(Y[Ta > SEQX[j]])

    if (Nb[j] < 11) {

      whitetest[j] <- NA
      p[j] <- NA
      p2[j] <- NA

    } else {
      whitetest[j] <- lmtest::bptest(stats::lm(Y[Ta > SEQX[j]] ~ Ta[Ta > SEQX[j]]))$p.value

      mod <- stats::lm(Y[Ta > SEQX[j]] ~ Ta[Ta> SEQX[j]])
      p[j] <- ifelse(stats::coefficients(mod)[2] > 0, NA, summary(mod)$coefficients[2,4])
      p2[j] <- ifelse(stats::coefficients(mod)[2] < 0, NA, summary(mod)$coefficients[2,4])


    }
  }

  if (length(SEQX[p2 < 0.01 & !is.na(p2)]) > 0 &
      length(SEQX[p < 0.01 & !is.na(p)]) == 0 &
      length(SEQX[whitetest < 0.05]) == 0) {

    stop("Tlc and Mtnz can not be estimated: provide values")

  } else {

    test_1 <- ifelse(length(stats::na.omit(SEQX[whitetest > 0.05])) >= 1, min(SEQX[whitetest > 0.05], na.rm = TRUE), NA)

    test_2 <- ifelse(length(stats::na.omit(SEQX[p < 0.01])) >= 1, SEQX[match(max(SEQX[p < 0.01],na.rm = TRUE) ,SEQX) - 1], NA)

    low_Tlc <- ifelse(length(stats::na.omit(c(test_1, test_2))) == 0, NA, max(test_1, test_2, na.rm = TRUE)) ### check if ifelse necessary

  }
  if (is.na(low_Tlc)) stop("Tlc and Mtnz can not be estimated: provide values")
  if (length(Y[Ta > low_Tlc]) < 10) warning("Mtnz computed on less than 10 points")

  Mtnz <- mean(Y[Ta > low_Tlc], na.rm = TRUE)
  out <- c(Mtnz,low_Tlc)
  names(out) <- c("Mtnz", "low_Tlc")
  out
}

#' Estimate Tlc
#'
#' This function estimate the lowest Tlc and the Mtnz
#'
#'@name estimate_Tlc_Mtnz
#'@param M A vector of methabolic measure.
#'@param Ta A vector of ambient temperature.
#'@param Mtnz A value of Mtnz.
#'@param fitting_options a list of fitting option to pass to jags.
#'@export
#'@examples
#'t <- estimate_Tlc_Mtnz(M = test_data2$VO2ms, Ta = test_data2$Ta, Mtnz = NULL,
#'                       fitting_options = list(parallel = FALSE))

estimate_Tlc_Mtnz <- function(M, Ta, Mtnz, fitting_options = list(ni = 50000,
                                                                  nt = 10,
                                                                  nb = 30000,
                                                                  nc = 3,
                                                                  parallel = TRUE)) { ## add parallel to fitting options
  .complete_args(estimate_Tlc_Mtnz)
  Y <- M
  Mtnz_input <- Mtnz


  out_1 <- find_low_Tlc_Mtnz(Y, Ta) ## run first step

  if(!is.null(Mtnz_input)) {
    Mtnz <- Mtnz_input
  } else {
    Mtnz <- out_1["Mtnz"]
  } ## extract Mtnz
  low_Tlc <- out_1["low_Tlc"] ## extract low Tlc


  Ym <- mean(Y, na.rm = TRUE) ## keep the mean
  Ta2 <- Ta[Ta < low_Tlc]
  Y2 <- Y[Ta < low_Tlc]

  inits <- list(
    tauy2= stats::runif(1,0.05,0.1),
    tauy1= stats::runif(1,0.03,0.05),
    p = stats::runif(2),
    Tbe= stats::runif(1,35,40),
    Tlc = stats::runif(1, low_Tlc, max(Ta)),
    Mr=stats::runif(1,Mtnz*0.8/Ym,Mtnz/Ym),
    TMR=stats::runif(1,0.7*Mtnz/(Ym),Mtnz*0.8/Ym))

  inits_hetero_list <- rep(list(inits), fitting_options[["nc"]])

  params_hetero <- c("Tlc")

  win_data <- list(Y = Y2/Ym,
                   NbObservations = length(Y2),
                   Ta = Ta2,
                   Mtnz = Mtnz/Ym,
                   Max = max(Ta))

  # MCMC settings
  path_to_model_1 <- system.file("extdata", "hetero1.txt",  package = "torpor")

  mod <- jagsUI::jags(data = win_data,
                      inits = inits_hetero_list,
                      parameters.to.save = params_hetero,
                      model.file = path_to_model_1,
                      n.chains = fitting_options[["nc"]],
                      n.thin = fitting_options[["nt"]],
                      n.iter = fitting_options[["ni"]],
                      n.burnin = fitting_options[["nb"]],
                      parallel = fitting_options[["parallel"]],
                      verbose = FALSE,
                      store.data = TRUE)

  ## extract the important output
  Tlc_estimated <- stats::median(mod$sims.list$Tlc) ## estimated Tlc
  Tlc_distribution <- mod$sims.list$Tlc ## distribution of Tlc

  if(is.null(Mtnz_input)) {
    Mtnz_estimated <- mean(Y[Ta >= Tlc_estimated], na.rm = TRUE) # estimated Mtnz (mean of the points)
    Mtnz_points <- stats::na.omit(Y[Ta >= Tlc_estimated])
  } else {
    Mtnz_estimated <- Mtnz_input
    Mtnz_points <- Mtnz_input
  }

  if(length(stats::na.omit(Y[Ta >= Tlc_estimated])) < 10 & is.null(Mtnz_input)) warning("Mtnz computed on less than 10 points")

  return(list(model_1 = mod,
              Tlc_estimated = Tlc_estimated,
              Tlc_distribution = Tlc_distribution,
              Mtnz_estimated = Mtnz_estimated,
              Mtnz_points = Mtnz_points,
              Ta2 = Ta2))
}

#' Estimate assignation
#'
#' This function estimate the assignation into torpor or euthermia. It is used
#' internally prior to estimate the parameters
#'
#'@name estimate_assignation
#'@export
#'@inheritParams estimate_Tlc_Mtnz
#'@param Mtnz Mtnz value if not estimated
#'@param Tlc Tlc value if not estimated
#'@examples
#'\dontrun{
#'t2_no_input <- estimate_assignation(Ta = test_data2$Ta, M = test_data2$VO2)
#'to_input <- estimate_assignation(Ta = test_data2$Ta, M = test_data2$VO2, Mtnz = 1.49, Tlc = 28.8)
#'}
estimate_assignation <- function(M, Ta,
                                 Mtnz = NULL, Tlc = NULL,
                                 fitting_options = list(ni = 50000,
                                                        nt = 10,
                                                        nb = 30000,
                                                        nc = 3,
                                                        parallel = TRUE)){


  .complete_args(estimate_assignation)

  ## check input
  Y <- M
  if (length(Y) != length(Ta)) stop("M and Ta have not the same length")

  data <- cbind(Y, Ta)[!is.na(Y) & !is.na(Ta), ] ## remove Nas

  Y <- data[, "Y"]
  Ta <- data[, "Ta"]
  Ym <- mean(Y, na.rm = T)

  ## run step 1 and 2 to get Mtnz and Tlc

  if (is.null(Mtnz) & (is.null(Tlc))) {
    message("Mtnz and Tlc are being estimated from the data")

    out_Mtnz_Tlc <- estimate_Tlc_Mtnz(Ta = Ta, M = Y, Mtnz = Mtnz, fitting_options = fitting_options)

  } else if (!(is.null(Mtnz)) & (is.null(Tlc))) {

    message("Tlc is being estimated from the data")

    out_Mtnz_Tlc <- estimate_Tlc_Mtnz(Ta = Ta, M = Y, Mtnz = Mtnz, fitting_options = fitting_options)

  } else if (is.null(Mtnz) & !(is.null(Tlc))) {

    message("Mtnz is being estimated from the data")

    out_Mtnz_Tlc <- list(model_1 = NULL,
                         Tlc_estimated = Tlc,
                         Tlc_distribution = Tlc,
                         Mtnz_estimated = mean(Y[Ta >= Tlc], na.rm = TRUE),
                         Mtnz_points = stats::na.omit(Y[Ta >= Tlc]),
                         Ta2 = stats::na.omit(Ta[Ta < Tlc]))

  } else { ## keep the structure of output of estimate_Tlc_Mtnz when there are provided as input
    message("test")
    out_Mtnz_Tlc <- list(model_1 = NULL,
                         Tlc_estimated = Tlc,
                         Tlc_distribution = Tlc,
                         Mtnz_estimated = Mtnz,
                         Mtnz_points = Mtnz,
                         Ta2 = Ta[Ta < Tlc])
  }

  Tlc <- out_Mtnz_Tlc$Tlc_estimated
  Mtnz <- out_Mtnz_Tlc$Mtnz_estimated

  ## initial values
  inits_2 <- list(
    tauy2=stats::runif(1,0.05,0.1),
    tauy1=stats::runif(1,0.03,0.05),
    p = stats::runif(3),
    Tbe = stats::runif(1,Tlc, 50),
    TMR = stats::runif(1,0.7*Mtnz/(Ym),Mtnz*0.8/Ym),
    Mr = stats::runif(1, Mtnz*0.8/Ym, Mtnz/Ym))


  inits_hetero_list_2 <- rep(list(inits_2), fitting_options[["nc"]])

  params_hetero_2 <- c("G", "Tt","intr","intc", "inte", "betat", "betac")
  Tlcq2.5 <- if(is.null(out_Mtnz_Tlc$model_1)) Tlc else out_Mtnz_Tlc$model_1$q2.5$Tlc

  win.data_2 <- list(Y = Y/Ym,
                     NbObservations = length(Y),
                     Ta = Ta,
                     Mtnz = Mtnz/Ym,
                     Tlc = Tlc,
                     Tlcq2.5= Tlcq2.5)

  path_to_model_2 <- system.file("extdata", "hetero2.txt",  package = "torpor")

  mod_assignation <- jagsUI::jags(data = win.data_2,
                                  inits = inits_hetero_list_2,
                                  parameters.to.save = params_hetero_2,
                                  model.file = path_to_model_2,
                                  n.chains = fitting_options[["nc"]],
                                  n.thin = fitting_options[["nt"]],
                                  n.iter = fitting_options[["ni"]],
                                  n.burnin = fitting_options[["nb"]],
                                  parallel = fitting_options[["parallel"]],
                                  verbose = FALSE,
                                  store.data = TRUE)

  out_assignation <- list(
    mod_assignation = mod_assignation,
    out_Mtnz_Tlc = out_Mtnz_Tlc,
    data = list(Y = Y,
                Ta = Ta,
                Ym = Ym))

  out_assignation$assignation <- .step_3_bis(out_assignation = out_assignation,
                                             fitting_options = fitting_options)

  out_assignation
}


#'step_3_bis
#'
#' This function estimate the lowest TLC and the Mtnz
#'
#'@name .step_3_bis
#'@param out_assignation output of estimatate_assignation.
#'@param fitting_options fitting options
.step_3_bis <- function(out_assignation, fitting_options){
  ## comes from hetero2
  Ta <- out_assignation$data$Ta
  Y <- out_assignation$data$Y

  nc <- fitting_options[["nc"]]
  nt <- fitting_options[["nt"]]
  ni <- fitting_options[["ni"]]
  nb <- fitting_options[["nb"]]
  ## extract model assignation
  mod_assignation <- out_assignation$mod_assignation

  G <- mod_assignation$mean$G ## extract the G


  nbsamples <- ((ni - nb)*nc)/nt ## look how to take them from mod2

  ## inside small fun
  expit <- function(x){1/(1 + exp(-x))} ##step not to be done if Tlc and Mtnz and given
  funabove <- function(x, mod){expit(stats::coefficients(mod)[1] + stats::coefficients(mod)[2]*x)} ##step not to be done if Tlc and Mtnz and given
  funbelow <- function(x, mod){-(funabove(x, mod)- 1)} ##step not to be done if Tlc and Mtnz and given

  ####
  if(length(round(G)[round(G)==1])!=0){

    for(i in 1:length(G)){
      #Only proceed to this step if there are at least one torpid (G==1) value.
      G[i] <- ifelse(Y[i] <= stats::median(tor_predict_fun(x = Ta[i],
                                                           Tt = mod_assignation$sims.list$Tt,
                                                           intr = mod_assignation$sims.list$intr,
                                                           intc = mod_assignation$sims.list$intc,
                                                           betat = mod_assignation$sims.list$betat,
                                                           betac = mod_assignation$sims.list$betac,
                                                           Ym = out_assignation$data$Ym)) & G[i] != 3, 1, G[i])

      G[i] <- ifelse(Y[i] >=  stats::median(eut_predict_fun(x = Ta[i],
                                                            inte = mod_assignation$sims.list$inte,
                                                            betat = mod_assignation$sims.list$betat,
                                                            Ym = out_assignation$data$Ym)) & G[i] != 3, 2 ,G[i])
    }
  }

  probStatus <- 1 - sqrt((G - round(G))^2)

  length_Tlc_dist <- length(out_assignation$out_Mtnz_Tlc$Tlc_distribution)

  if(length_Tlc_dist > 1) {
    mod_glm <- stats::glm(stats::pnorm(q = out_assignation$out_Mtnz_Tlc$Tlc_distribution, ## changed for Tlc distribution
                                       mean = out_assignation$out_Mtnz_Tlc$Tlc_estimated,
                                       sd = ifelse(length_Tlc_dist > 1, stats::sd(out_assignation$out_Mtnz_Tlc$Tlc_distribution), 0)) ~
                            out_assignation$out_Mtnz_Tlc$Tlc_distribution,
                          family = "binomial") ##step not to be done if Tlc and Mtnz and given


    probStatus <- probStatus*ifelse(G == 3, funabove(Ta, mod_glm), funbelow(Ta, mod_glm))  ##step not to be done if Tlc and Mtnz and given
  }
  pstatus <- rep(NA,length(Ta))

  for(i in 1:length(Ta)) {

    pstatus[i] <- as.numeric(stats::binom.test(round(nbsamples*probStatus[i]),
                                               nbsamples,
                                               alternative = "greater", p = 0.5)$p.value)

  }

  G[pstatus>0.01] <- 0
  G <- round(G)



  list(G = G,
       pstatus = pstatus)

}


#' Mixture model aimed at assigning metabolic rate measurements (M) to torpor and euthermia
#'
#'The function [tor_fit()] considers the relation between metabolic rate (M) and
#'ambient temperature (Ta) assumed by the Scholander-Irving model and its later extensions.
#'
#'Resting M measured within the thermoneutral zone (TNZ) is independent of Ta.
#'This rate is hereafter referred to as Mtnz, although it would correspond to
#'the basal metabolic rate (BMR) provided that the specific criteria for the
#'BMR are met (see Fasel et al. in prep.). Below the lower critical temperature
#'of TNZ (Tlc), M of euthermic animals increases linearly with decreasing Ta.
#'M of torpid animals increases linearly with decreasing Ta to maintain a minimal
#'body temperature below some threshold ambient temperature (Tt). This state is
#'usually referred to as "regulated torpor". Between Tt and Tlc, M of torpid
#'animals follows an exponential curve. In this Ta range, torpor is referred to
#'as "conforming torpor".
#'
#'@name tor_fit
#'@inheritParams estimate_assignation
#'@return A list of class tor_obj
#'@export
#'@examples
#'\dontrun{
#'test_mod <- tor_fit(Ta = test_data$Ta, M = test_data$VO2,
#'                    fitting_options = list(parallel = TRUE))
#'test_mod2 <- tor_fit(Ta = test_data2$Ta, M = test_data2$VO2ms, Mtnz = 1.8, Tlc = 29.2,
#'                    fitting_options = list(parallel = TRUE))
#'test_mod3 <- tor_fit(Ta = test_data2$Ta, M = test_data2$VO2ms, Mtnz = 1.8,
#'                    fitting_options = list(parallel = TRUE))
#'
#'test_mod4 <- tor_fit(Ta = test_data2$Ta, M = test_data2$VO2ms, Tlc = 29.2,
#'                    fitting_options = list(parallel = TRUE))
#'}
tor_fit <- function(Ta, M,
                    Mtnz = NULL,
                    Tlc = NULL,
                    fitting_options = list(ni = 50000,
                                           nt = 10,
                                           nb = 30000,
                                           nc = 3,
                                           parallel = TRUE)) {


  if (length(Ta) != length(M)) stop("M and Ta do not have the same length")
  ## run step 1-3
  .complete_args(tor_fit)


  out_assignation <- estimate_assignation(M = M, Ta = Ta, Mtnz = Mtnz, Tlc = Tlc, fitting_options)

  G <- out_assignation$assignation$G
  Y <- out_assignation$data$Y
  Ta <- out_assignation$data$Ta
  Tlc <- out_assignation$out_Mtnz_Tlc$Tlc_estimated
  Mtnz <- out_assignation$out_Mtnz_Tlc$Mtnz_estimated
  Ym <- out_assignation$data$Ym

  Tlcq2.5 <- if(is.null(out_assignation$out_Mtnz_Tlc$model_1)) Tlc else out_assignation$out_Mtnz_Tlc$model_1$q2.5$Tlc ##

  win.data_3 <- list(Y = Y[G != 0 & G != 3] / Ym,
                     NbObservations = length(Y[G != 0 & G != 3]),
                     Ta = Ta[G != 0 & G != 3],
                     Mtnz = Mtnz / Ym,
                     Tlc = Tlc,
                     G = G[G != 0 & G != 3],
                     Tlcq2.5= Tlcq2.5)

  path_to_model_3 <- system.file("extdata", "hetero3.txt",  package = "torpor")

  inits_3 <-  list(
    tauy2=stats::runif(1,0.05,0.1),
    tauy1=stats::runif(1,0.03,0.05),
    Tbe = stats::runif(1,Tlc, 50),
    TMR = stats::runif(1,0.7*Mtnz/(Ym),Mtnz*0.8/Ym),
    Mr = stats::runif(1, Mtnz*0.8/Ym, Mtnz/Ym))

  inits_hetero_list_3 <- rep(list(inits_3), fitting_options[["nc"]])
  params <- params_hetero_2 <- c("tau","inte","intc","intr","betat",
                                 "betac","Tt","TMR","Mr", "Tbt", "Tbe")

  out_4 <- jagsUI::jags(data = win.data_3,
                        inits =  inits_hetero_list_3,
                        parameters.to.save = params,
                        model.file = path_to_model_3,
                        n.chains = fitting_options[["nc"]],
                        n.thin = fitting_options[["nt"]],
                        n.iter = fitting_options[["ni"]],
                        n.burnin = fitting_options[["nb"]],
                        parallel = fitting_options[["parallel"]],
                        verbose = FALSE,
                        store.data = TRUE)

  ### check convergence.... if Rhat > 1.1 warrning sur parameters.
  for(i in seq_along(out_4$Rhat)){
    if (any(out_4$Rhat[i][[1]] > 1.1, na.rm = TRUE)) warning(paste("Rhat > 1.1 on parameter:", names(out_4$Rhat[i])))
  }

  out_assignation$mod_parameter <- out_4

  class(out_assignation) <- c("tor_obj", "list")
  print(out_assignation)
  return(invisible(out_assignation))

}
