#' Fit Torpor
#'
#'wrapper around the Jags function -
#'@name fit_torpor
#'@param MR a vector of Metabolic rate
#'@param Ta a vector of Temperature (same length than Ta)
#'@param BMR value for the focal specie
#'@param TLC value for the focal specie
#'@param Model path to model_file.txt
#'@param fitting_options a list with specification for jags parameters
#'@return a fitted jags object
#'@import rjags
#'@import jagsUI
#'@export
#'
fit_torpor <- function(MR,
                       Ta,
                       BMR,
                       TLC,
                       Model = NULL,
                       fitting_options = list(ni = 500000,
                                              nt = 10,
                                              nb = 300000,
                                              nc = 3)){

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
                     parallel=T)

return(out)
}

