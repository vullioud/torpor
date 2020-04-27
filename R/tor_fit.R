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
#' This function estimate the lowest TLC and
#'
#'@name complete_args
#'@aliases complete_args
#'@param Y A vector of methabolic measure.
#'@param Ta A vector of ambient temperature.
#'@examples
#'step_1(test_data$VO2, test_data$Ta)

step_1 <- function(MR,Ta){

  SEQX <- seq(from = max(Ta,na.rm = TRUE), to = min(Ta,na.rm = TRUE), length.out = 100)
  Nb <- rep(NA, 100)
  whitetest <- rep(NA, 100)
  p <- rep(NA, 100)

  #Define heterosc and slope
  for (j in 1:100) {

      Nb[j] <- length(MR[Ta > SEQX[j]])
    if (Nb[j] < 11) {
      whitetest[j] <- NA
    p[j] <- NA
    } else {
      whitetest[j] <- lmtest::bptest(lm(MR[Ta > SEQX[j]] ~ Ta[Ta > SEQX[j]]))$p.value

      mod <- lm(MR[Ta> SEQX[j]]~Ta[Ta> SEQX[j]])
      p[j] <- ifelse(coefficients(mod)[2] > 0,NA,summary(mod)$coefficients[2,4])}
  }

  #Define "low TLC" and first BMR

  test_1 <- ifelse(length(na.omit(SEQX[whitetest > 0.05])) >= 1, min(SEQX[whitetest > 0.05], na.rm = TRUE), NA)
  test_2 <- ifelse(length(na.omit(SEQX[p < 0.01])) >= 1, SEQX[match(max(SEQX[p < 0.01],na.rm = TRUE) ,SEQX) - 1], NA)

  heterosc <- ifelse(length(na.omit(c(test_1, test_2))) < 1, NA, max(test_1, test_2, na.rm = T))

  if (is.na(heterosc)) stop("TLC and BMR not estimable, please provide them by hand")

  bmr <- mean(MR[Ta > heterosc], na.rm = TRUE)
out <- c(bmr,heterosc)
names(out) <- c("bmr", "low_tlc")
return(out)
}
