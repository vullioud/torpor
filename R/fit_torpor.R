#' Fit Torpor
#'
#'The function \code{fit_torpor} fit binomial mixture model using Bayesian inference.
#'This function uses Rjags in the background and enable users to specify some
#'sampling parameters and, if wanted, change the structure of the model directly.
#'
#'The function considers the assumed relation between MR and Ta (Speakman & Thomas 2003).
#'In the hypothermic state (torpor) and above some threshold Ta (Tmin),
#'MR follows an exponential curve reflecting the Arrhenius rate enhancing effect
#'of temperature on chemical reactions, whereas below Tmin, it increases linearly
#'with decreasing Ta to maintain a minimal Tb in torpor.
#'In the euthermic state, MR solely increases linearly with decreasing Ta.
#'
#'@name fit_torpor
#'@param MR a vector of Metabolic rate
#'@param Ta a vector of ambient Temperature (same length as MR)
#'@param BMR BMR value for the focal specie
#'@param TLC TLC value for the focal specie
#'@param Model path to model_file.txt if the users want to use her/his own model structure
#'@param fitting_options a list with specification for sampling parameters.
#'The follwing parameters can be speficied:
#'*ni = number of interation
#'*nt = number of XXX
#'*nb = number of burns
#'*nc = number of chain
#'
#'@return a fitted jags object
#'@import rjags
#'@import jagsUI
#'@export
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
                                              nc = 3)){

  CompleteArgs(fit_torpor)
## find the model if not given by the user
 if(is.null(Model)){
   path_to_model <- system.file("extdata", "hetero.txt",  package = "toRpoR")
 } else {
   path_to_model <- paste(Model)
  }

## Returned values
params_hetero <- c("tauy",
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
                   "TLC")

## get the values for the models
Y <- as.numeric(as.character(MR))
Ta <- as.numeric(as.character(Ta))

## check
if(length(Y) != length(Ta)) {
  stop("Ta and MR not the same length")
}

## remove NAs
da <- cbind(Y,Ta)[!is.na(Y)&!is.na(Ta)&Ta<(TLC-2),]
Y <-as.numeric(da[,1])
Ta <-as.numeric(da[,2])


## shuffle the data
set.seed(666)
ord <- order(stats::rnorm(length(Y),0,1))
Y <- Y[ord]
Ta <- Ta[ord]

# create data list
win.data <- list(Y = Y,
                 NbObservations=length(Y),
                 Ta=Ta,
                 TLC=TLC,
                 BMR=BMR)

# create initial values
inits<- list(
    tauy=stats::runif(1),
    p=stats::runif(2),
    beta1=stats::runif(1,-0.1,0),
    beta2=stats::runif(1,0.05,0.1),
    TMR=stats::runif(1,0.1,0.18),
    tlc=stats::runif(1,TLC,40)
  )
inits_hetero_list <- rep(list(inits), fitting_options[["nc"]])

out <- jagsUI ::jags(data = win.data,
                     inits = inits_hetero_list,
                     parameters.to.save = params_hetero,
                     model.file = path_to_model,
                     n.chains = fitting_options[["nc"]],
                     n.thin = fitting_options[["nt"]],
                     n.iter = fitting_options[["ni"]],
                     n.burnin = fitting_options[["nb"]],
                     parallel=T, verbose = F)

return(out)
}


#' CompleteArgs
#'
#'From Alex Courtiol Isorix package
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
