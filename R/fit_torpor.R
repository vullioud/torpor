#' Fit Torpor
#'
#'The function fit_torpor() fits a binomial mixture model using Bayesian
#'inference. It uses Rjags in the background and enables users to specify
#'some - but not all - sampling parameters. The structure of the model can also
#'be changed. User who want more flexibility are encouraged to use Rjags
#'directly. The function considers the assumed relation between MR and Ta
#'(Speakman & Thomas 2003). In the hypothermic state (torpor) and above some
#'threshold Ta (Tmin), MR follows an exponential curve reflecting the Arrhenius
#'rate enhancing effect of temperature on chemical reactions, whereas below
#'Tmin, it increases linearly with decreasing Ta to maintain a minimal Tb in
#'torpor. In the euthermic state, MR solely increases linearly with decreasing Ta.
#'
#'This model should be applied only with sufficient sample size and if
#'evidences suggest that individuals under study will conform to the previously
#'described pattern while in torpor.
#'
#'@name fit_torpor
#'@aliases fit_torpor
#'@param MR a vector of Metabolic rates
#'@param Ta a vector of ambient Temperatures (same length as MR)
#'@param BMR BMR value for the focal specie
#'@param TLC TLC value for the focal specie
#'@param Model path to model_file.txt if a different model is used
#'@param fitting_options a list specifying sampling parameters. The follwing parameters can be speficied:
#' * ni = number of itterations: default ni = 500000
#' * nt = thin rate: default nt = 10
#' * nb = number of burns: default nb = 300000
#' * nc = number of chains: default nc = 3
#'
#'@return fitted model. A list of class "JagsUI".
#'@export
#'@import rjags
#'@import jagsUI
#'@importFrom stats runif
#'@examples
#'rm(list = ls())
#'data(test_data3)
#'test2 <- fit_torpor(MR = test_data3[,2],
#'Ta = test_data3[, 1],
#'BMR = 1.005,
#'TLC = 29,
#'fitting_options = list(nc = 2))

fit_torpor <- function(MR,
                       Ta,
                       BMR,
                       TLC,
                       Model = NULL,
                       fitting_options = list(ni = 500000,
                                              nt = 10,
                                              nb = 300000,
                                              nc = 3)) {
  complete_args(fit_torpor)

  ## check input
  if (length(MR) != length(Ta)) {
    stop("Ta and MR have not the same length")
  }

  ## find the model if not given by the user
  if (is.null(Model)) {
    path_to_model <- system.file("extdata", "hetero.txt",  package = "toRpoR")
  } else {
    path_to_model <- paste(Model)
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
                   NbObservations=length(Y),
                   Ta=Ta,
                   TLC = TLC,
                   BMR = BMR/Ym,
                   Ym = Ym) ## included to be able to back-transform BMR and Y for predictions.


  # create initial values
  inits <- list(
      tau=runif(1,0.1,1),
      p=runif(2),
      betat=runif(1,-0.1,0),
      tlc =runif(1,TLC,TLC+1),
      Tt=runif(1),
      TMR=runif(1,0,0.1),
      diff=runif(1,1,2))

  # create initial values
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
#'complete_args
#'
#'this function is used internally for fit_torpor
#'From Alex Courtiol
#'
#'@name complete_args
#'@aliases complete_args
#'@param fn a function

complete_args <- function(fn) {
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
