#' tor_fit
#'
#'[tor_fit()] fits a binomial mixture model using Bayesian
#'inference.  The function considers the assumed relation between metabolic rate (MR) and ambient temperature (Ta)
#'(Speakman & Thomas 2003). In the hypothermic state (torpor) and above some
#'threshold Ta (Tmin), MR follows an exponential curve reflecting the Arrhenius
#'rate enhancing effect of temperature on chemical reactions, whereas below
#'Tmin, it increases linearly with decreasing Ta to maintain a minimal body temperature (Tb) in
#'torpor. In the euthermic state, MR solely increases linearly with decreasing Ta.
#'The function uses Rjags in the background and enables users to specify
#'some - but not all - sampling parameters. The structure of the model can also
#'be changed. Users who want more flexibility are encouraged to use Rjags
#'directly.
#'More information about the model can be found in the \code{vignette("model_description", package = "torpor")}.
#'
#'This model should be applied with sufficient sample size and if
#'evidence suggest that individuals under study will conform to the previously
#'described pattern while in torpor.
#'
#'@name tor_fit
#'@aliases tor_fit
#'@family fit
#'@param MR a vector of Metabolic rates
#'@param Ta a vector of ambient Temperatures (same length as MR)
#'@param BMR BMR Basal Metabolic Rate for the focal specie
#'@param TLC lower critical temperature for the focal specie
#'@param model path to model_file.txt if a different model is used
#'@param fitting_options a list specifying sampling parameters. The follwing parameters can be speficied:
#' * ni = number of itterations: default ni = 500000
#' * nt = thin rate: default nt = 10
#' * nb = number of burns: default nb = 300000
#' * nc = number of chains: default nc = 3
#'@return fitted model. A list of class "JagsUI".
#'@export
#'@importFrom stats runif
#'@examples
#'rm(list = ls())
#'data(test_data3)
#'test2 <- tor_fit(MR = test_data3[,2],
#'Ta = test_data3[, 1],
#'BMR = 1.005,
#'TLC = 29,
#'model = NULL,
#'fitting_options = list(ni = 5000, nb = 3000, nc = 2))

tor_fit <- function(MR,
                    Ta,
                    BMR,
                    TLC,
                    model = NULL,
                    fitting_options = list(ni = 500000,
                                           nt = 10,
                                           nb = 300000,
                                           nc = 3)) {
  .complete_args(tor_fit)

  ## check input
  if (length(MR) != length(Ta)) {
    stop("Ta and MR have not the same length")
  }

  ## find the model if not given by the user
  if (is.null(model)) {
    path_to_model <- system.file("extdata", "hetero.txt",  package = "torpor")
  } else {
    path_to_model <- paste(model)
  }

  ## Returned values
  params_hetero <- c("tauy",
                     "G",
                     "p",
                     "inte",
                     "intc",
                     "intr",
                     "betat",
                     "betac",
                     "Tt",
                     "TMR",
                     "tlc",
                     "diff",
                     "TLC",
                     "BMR",
                     "Ym")

  ## get the values for the models /reorder the data and remove NA
  da <- cbind(MR, Ta)[!is.na(MR) & !is.na(Ta) & Ta < (TLC), ]
  Y <- da[, "MR"]
  Ta <- da[, "Ta"]

  ## get the mean to standardised Y and BMR
  Ym <- mean(Y, na.rm = TRUE)

  # create data list
  win.data <- list(Y = Y/Ym,
                   NbObservations = length(Y),
                   Ta = Ta,
                   TLC = TLC,
                   BMR = BMR/Ym,
                   Ym = Ym) ## included to be able to back-transform BMR and Y for predictions.


  # create initial values
  inits <- list(
      tau = runif(1,0.1,1),
      p = runif(2),
      betat = runif(1,-0.1,0),
      tlc = runif(1,TLC,TLC + 1),
      Tt = runif(1),
      TMR = runif(1,0,0.1),
      diff = runif(1,1,2))

  inits_hetero_list <- rep(list(inits), fitting_options[["nc"]])


  out <- jagsUI::jags(
    data = win.data,
    inits = inits_hetero_list,
    parameters.to.save = params_hetero,
    model.file = path_to_model,
    n.chains = fitting_options[["nc"]],
    n.thin = fitting_options[["nt"]],
    n.iter = fitting_options[["ni"]],
    n.burnin = fitting_options[["nb"]],
    parallel = TRUE,
    verbose = FALSE,
    store.data = TRUE
  )

  return(out)
}

##############################################################################

#'complete_args // internal
#'
#'this function is used internally for tor_fit
#'From Alex Courtiol - informed. He does not want to be listed as an author -
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

#'step_1
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_1
#'@param Y A vector of methabolic measure.
#'@param Ta A vector of ambient temperature.
#'@export
#'@examples
#'step_1(test_data$VO2, test_data$Ta)

step_1 <- function(Y,Ta){

  SEQX <- seq(from = max(Ta,na.rm = TRUE), to = min(Ta,na.rm = TRUE), length.out = 100)
  Nb <- rep(NA, 100)
  whitetest <- rep(NA, 100)
  p <- rep(NA, 100)

  #Define heterosc and slope
  for (j in 1:100) {

      Nb[j] <- length(Y[Ta > SEQX[j]])
    if (Nb[j] < 11) {
      whitetest[j] <- NA
    p[j] <- NA
    } else {
      whitetest[j] <- lmtest::bptest(stats::lm(Y[Ta > SEQX[j]] ~ Ta[Ta > SEQX[j]]))$p.value

      mod <- stats::lm(Y[Ta> SEQX[j]]~Ta[Ta> SEQX[j]])
      p[j] <- ifelse(stats::coefficients(mod)[2] > 0,NA, summary(mod)$coefficients[2,4])}
  }

  #Define "low TLC" and first BMR

  test_1 <- ifelse(length(stats::na.omit(SEQX[whitetest > 0.05])) >= 1, min(SEQX[whitetest > 0.05], na.rm = TRUE), NA)
  test_2 <- ifelse(length(stats::na.omit(SEQX[p < 0.01])) >= 1, SEQX[match(max(SEQX[p < 0.01],na.rm = TRUE) ,SEQX) - 1], NA)

  heterosc <- ifelse(length(stats::na.omit(c(test_1, test_2))) < 1, NA, max(test_1, test_2, na.rm = TRUE))

  if (is.na(heterosc)) stop("TLC and BMR not estimable, please provide them by hand")

  bmr <- mean(Y[Ta > heterosc], na.rm = TRUE)
  out <- c(bmr,heterosc)
names(out) <- c("bmr", "low_tlc")
return(out)
}

#'step_1_and_2
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_1_and_2
#'@param Y A vector of methabolic measure.
#'@param Ta A vector of ambient temperature.
#'@inheritParams tor_fit
#'@export
step_1_and_2 <- function(Y, Ta, fitting_options = list(ni = 500000,
                                                       nt = 10,
                                                       nb = 300000,
                                                       nc = 3)) {

  out_1 <- step_1(Y, Ta) ## run first step

  bmr <- out_1["bmr"] ## extract bmr
  heterosc <- out_1["low_tlc"] ## extract low tlc

  Ym <- mean(Y, na.rm = TRUE) ## keep the mean
  Ta2 <- Ta[Ta < heterosc]
  Y2 <- Y[Ta < heterosc]

  inits <- list(
    tauy = runif(1,0.05,1),
    coef = runif(1,1,1.1),
    p = runif(2),
    betat = runif(1,-0.1,0),
    tlc = runif(1,heterosc,max(Ta)),
    Tt = runif(1),
    TMR = runif(1,0,bmr*0.8/Ym),
    intc1 = runif(1,0,bmr*0.8/Ym))

  inits_hetero_list <- rep(list(inits), fitting_options[["nc"]])

  params_hetero <- c("tlc","betat","inte")

  win_data <- list(Y = Y2/Ym,
                   NbObservations = length(Y2),
                   Ta = Ta2,
                   heterosc = heterosc,
                   BMR = bmr/Ym,
                   Max = max(Ta),
                   Ym = Ym)

  # MCMC settings


  path_to_model_1 <- system.file("extdata", "hetero_1.txt",  package = "torpor")

  out <- jagsUI::jags(data = win_data,
                        inits = inits_hetero_list,
                        parameters.to.save = params_hetero,
                        model.file = path_to_model_1,
                        n.chains = fitting_options[["nc"]],
                        n.thin = fitting_options[["nt"]],
                        n.iter = fitting_options[["ni"]],
                        n.burnin = fitting_options[["nb"]],
                        parallel = TRUE,
                        verbose = FALSE,
                        store.data = TRUE)

return(out)
}

#'step_1_and_2
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_3
#'@export
#'@param bmr bmr value if not estimated
#'@param tlc tlc value if not estimated
#'@inheritParams step_1_and_2
#'@examples
#'\dontrun{
#'step_3(Ta = test_data$Ta, Y = test_data$VO2)
#'}
step_3 <- function(Y, Ta, bmr = NULL, tlc = NULL, fitting_options = list(ni = 500000,
                                                                         nt = 10,
                                                                         nb = 300000,
                                                                         nc = 3)){

  .complete_args(step_3)

  ## check input
  if (length(Y) != length(Ta)) {
    stop("Ta and MR have not the same length")
  }

  data <- cbind(Y, Ta)[!is.na(Y) & !is.na(Ta), ] ## remove Nas

  Y <- data[, "Y"]
  Ta <- data[, "Ta"]

  ## run step 1 and 2 to get bmr and tlc

  if (is.null(bmr) | is.null(tlc)) {
    message("BMR and TLC are being
            estimated from the data")
    out_2 <- step_1_and_2(Ta = Ta, Y = Y)
    bmr <- mean(Y[Ta > out_2$mean$tlc])
    tlc <- out_2$mean$tlc
  }

  Ym <- mean(Y, na.rm = T)

  ## initial values
  inits_2 <- list(
    tauy = runif(1,0.05,1),
    coef = runif(1,1,1.1),
    p = runif(3),
    betat = runif(1,-0.1,0),
    intc1 = runif(1,0, bmr*0.8/Ym),
    Tt = runif(1),
    TMR = runif(1,0, bmr*0.8/Ym))

  inits_hetero_list_2 <- rep(list(inits_2), fitting_options[["nc"]])

  params_hetero_2 <- c("tau",
                       "G",
                       "p",
                       "inte",
                       "intc",
                       "intr",
                       "betat",
                       "betac",
                       "Tt",
                       "TMR",
                       "tlc")

  win.data_2 <- list(Y = Y/Ym,
                     NbObservations = length(Y),
                     Ta = Ta,
                     BMR = bmr/Ym,
                     tlc = tlc,
                     Ym = Ym)

  path_to_model_2 <- system.file("extdata", "hetero_2.txt",  package = "torpor")

  out2 <- jagsUI::jags(data = win.data_2,
                       inits = inits_hetero_list_2,
                       parameters.to.save = params_hetero_2,
                       model.file = path_to_model_2,
                       n.chains = fitting_options[["nc"]],
                       n.thin = fitting_options[["nt"]],
                       n.iter = fitting_options[["ni"]],
                       n.burnin = fitting_options[["nb"]],
                       parallel = TRUE,
                       verbose = FALSE,
                       store.data = TRUE)

  return(out2)

}


#'step_3_bis
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_3_bis
#'@param mod A fitted model
#'@export
step_3_bis <- function(mod){

  betat <- mod$sims.list$betat
  betac <- mod$sims.list$betac
  inte <- mod$sims.list$inte
  intc <- mod$sims.list$intc
  intr <- mod$sims.list$intr
  Tt <- mod$sims.list$Tt

  Ym <- mod$data$Ym
  Ta <- mod$data$Ta
  Y <- mod$data$Y
  tlc <- mod$data$tlc

  G <- mod$mean$G

  G[sqrt((G - round(G))^2) > 0.496] <- 0
  G <- round(G)

  if (length(G[G == 1]) != 0 ) {

    for (i in 1:length(Y <= tlc)) {

      #Only proceed to this step if there are at least one torpid (G==1) value.
      G[i] <- ifelse(Y[i] <= stats::median(tor_predict_fun(Ta[i], Tt, intr, intc, betat, betac, Ym)) & G[i] != 3, 1, G[i])
      G[i] <- ifelse(Y[i] >=  stats::median(eut_predict_fun(Ta[i], inte, betat, Ym)) & G[i] != 3, 2 ,G[i])
    }

  }
  return(G)
}

#'step_4
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_4
#'@inheritParams step_3
#'@export
#'@examples
#'\dontrun{
#'test_mod <- step_4(Ta = test_data3$Ta, Y = test_data3$VO2)
#'}
step_4 <- function(Ta, Y,  bmr = NULL, tlc = NULL, fitting_options = list(ni = 500000,
                                                                          nt = 10,
                                                                          nb = 300000,
                                                                          nc = 3)) {


  ## run step 1-3
  out_2 <- step_3(Ta = Ta, Y = Y,  bmr = bmr, tlc = tlc, fitting_options = fitting_options)

  ## manual
  tt <- step_3_bis(out_2)

  Y <- out_2$data$Y
  Ta <- out_2$data$Ta
  Ym <- out_2$data$Ym
  tlc <- out_2$data$tlc
  bmr <- out_2$data$BMR

  G <- tt

  if (length(Ta) != length(Y) | length(Y) != length(G)) stop("Y, Ta and G don´t have the same length")

  win.data_3 <- list(Y = Y[G != 0 & G != 3] / Ym,
                     NbObservations = length(Y[G != 0 & G != 3]),
                     Ta = Ta[G != 0 & G != 3],
                     BMR = bmr / Ym,
                     tlc = tlc,
                     G = G[G != 0 & G != 3],
                     Ym = Ym)


  path_to_model_3 <- system.file("extdata", "hetero_3.txt",  package = "torpor")

  inits_3 <- list(
    tauy = runif(1,0.05,1),
    coef = runif(1,1,1.1),
    p = runif(3),
    betat = runif(1,-0.1,0),
    intc1 = runif(1,0, bmr*0.8/Ym),
    Tt = runif(1),
    TMR = runif(1,0, bmr*0.8/Ym))

  inits_hetero_list_3 <- rep(list(inits_3), fitting_options[["nc"]])

  out3 <- jagsUI::jags(data = win.data_3,
                       inits =  inits_hetero_list_3,
                       parameters.to.save = out_2$parameters,
                       model.file = path_to_model_3,
                       n.chains = fitting_options[["nc"]],
                       n.thin = fitting_options[["nt"]],
                       n.iter = fitting_options[["ni"]],
                       n.burnin = fitting_options[["nb"]],
                       parallel = TRUE,
                       verbose = FALSE,
                       store.data = TRUE)


  return(out3)
}
