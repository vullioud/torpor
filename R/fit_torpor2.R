#' Fit Torpor 2
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
#'@name fit_torpor2
#'@aliases fit_torpor2
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
#'@return jagsUI object
#'@export
#'@import rjags
#'@import jagsUI
#'@importFrom stats runif
#'@examples
#'data(test_data)
#'test <- fit_torpor2(MR = test_data[,2],
#'Ta = test_data[, 1],
#'BMR = 98,
#'TLC = 28.88,
#'Model = NULL,
#'fitting_options = list(nc = 3))

fit_torpor2 <- function(MR,
                       Ta,
                       BMR,
                       TLC,
                       Model = NULL,
                       fitting_options = list(ni = 50000,
                                              nt = 2,
                                              nb = 30000,
                                              nc = 3)) {
  CompleteArgs(fit_torpor2)

  ## check input
  if (length(MR) != length(Ta)) {
    stop("Ta and MR have not the same length")
  }

  ## find the model if not given by the user
  if (is.null(Model)) {
    path_to_model <- system.file("extdata", "hetero2.txt",  package = "toRpoR")
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
  set.seed(666)
  da <- cbind(MR, Ta)[!is.na(MR) & !is.na(Ta) & Ta < (TLC - 2), ]
  da <- da[sample(nrow(da)), ]
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

