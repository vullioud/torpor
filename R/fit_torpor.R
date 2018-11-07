#' Fit Torpor
#'
#'The function \code{fit_torpor} fit binomial mixture model using Bayesian inference.
#'It uses Rjags in the background and enables users to specify some - but not all -
#'sampling parameters. The structure of the model can also be changed. User who want more
#'flexibility are encouraged to use Rjags directly.
#'The function considers the assumed relation between MR and Ta (Speakman & Thomas 2003).
#'In the hypothermic state (torpor) and above some threshold Ta (Tmin),
#'MR follows an exponential curve reflecting the Arrhenius rate enhancing effect
#'of temperature on chemical reactions, whereas below Tmin, it increases linearly
#'with decreasing Ta to maintain a minimal Tb in torpor.
#'In the euthermic state, MR solely increases linearly with decreasing Ta.
#'CAUTION: This model should be applied only if enough evidence is available suggesting
#'that the individuals under study will conform to the previously described pattern while in torpor.
#'@name fit_torpor
#'@param MR a vector of Metabolic rates
#'@param Ta a vector of ambient Temperatures (same length as MR)
#'@param BMR BMR value for the focal specie
#'@param TLC TLC value for the focal specie
#'@param Model path to model_file.txt if a different model is used
#'@param fitting_options a list specifying sampling parameters. The follwing parameters can be speficied:
#' * ni = number of itterations
#' * nt = thin rate
#' * nb = number of burns
#' * nc = number of chains
#'@return a fitted jags object
#'@export
#'@import rjags
#'@import jagsUI
#'@examples
#'data(test_data)
#'fit_torpor(MR = test_data[,2],
#'Ta = test_data[, 1],
#'BMR = 98,
#'TLC = 28.88,
#'Model = NULL,
#'fitting_options = list(nc = 1))

fit_torpor <- function(MR,
                       Ta,
                       BMR,
                       TLC,
                       Model = NULL,
                       fitting_options = list(ni = 500000,
                                              nt = 10,
                                              nb = 300000,
                                              nc = 3)) {
  CompleteArgs(fit_torpor)

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
  params_hetero <- c(
    "tauy",
    "G",
    "p",
    "int1",
    "int2",
    "int3",
    "beta1",
    "beta2",
    "Tmin",
    "tlc",
    "BMR",
    "TMR",
    "TLC"
  )

  ## get the values for the models /reorder the data and remove NA
  set.seed(666)
  da <- cbind(MR, Ta)[!is.na(MR) & !is.na(Ta) & Ta < (TLC - 2), ]
  da <- da[sample(nrow(da)), ]
  Y <- da[, "MR"]
  Ta <- da[, "Ta"]

  # create data list
  win.data <- list(
    Y = Y,
    NbObservations = length(Y),
    Ta = Ta,
    TLC = TLC,
    BMR = BMR
  )

  # create initial values
  inits <- list(
    tauy = stats::runif(1),
    p = stats::runif(2),
    beta1 = stats::runif(1, -0.1, 0),
    beta2 = stats::runif(1, 0.05, 0.1),
    TMR = stats::runif(1, 0.1, 0.18),
    tlc = stats::runif(1, TLC, 40)
  )
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
    parallel = T,
    verbose = F
  )

  return(out)
}


#' CompleteArgs
#'
#'borrowed from Alex Courtiol Isorix package
#'@name CompleteArgs
#'@param fn a function
CompleteArgs <- function(fn) {
  ## This function should not be called by the user but is itself called by other functions.
  ## It keeps the default list elements when
  ## a new list with fewer elements is provided
  env <- parent.frame()
  args <- formals(fn)
  for (arg.name in names(args)) {
    if (is.call(arg <- args[[arg.name]])) {
      if (arg[1] == "list()") {
        arg.input <- mget(names(args), envir = env)[[arg.name]]
        arg.full  <- eval(formals(fn)[[arg.name]])
        arg.full.updated <- utils::modifyList(arg.full, arg.input)
        assign(arg.name, arg.full.updated, envir = env)
      }
    }
  }
  return(NULL)
}
