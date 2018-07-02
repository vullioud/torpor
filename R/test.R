#' Fit models
#'
#'wrapper around the Jags function -
#'matic rounding
#'@name fit_model
#'@param data a data frame with 2 col. 1: Temperature, & 2: Metabolic value
#'@param BMR value for the focal specie
#'@param TLC value for the focal specie
#'@param Model path to model file
#'@param fitting_options a list with specification for jags parameters
#'@return a fitted jags object
#'@import rjags
#'@import dplyr

fit_torpor <- function(data, BMR, TLC, Model = NULL, fitting_options = list(ni = 500000,
                                                                            nt = 10,
                                                                            nb = 300000,
                                                                            nc = 3)){
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
                     "TMR")

## get the values for the models
Y <- as.numeric(as.character(data[,1]))
Ta <- as.numeric(as.character(data[,2]))
## remove NAs
da <- cbind(Y,Ta)[!is.na(Y)&!is.na(Ta)&Ta<(TLC-2),]
Y <-as.numeric(da[,1])
Ta <-as.numeric(da[,2])


## shuffle the data
set.seed(666)
ord <- order(rnorm(length(Y),0,1))
Y <- Y[ord]
Ta <- Ta[ord]

#
win.data <- list(Y = Y,
                 NbObservations=length(Y),
                 Ta=Ta,
                 TLC=Tlc,
                 BMR=BMR)

out <- jagsUI ::jags(win.data,
                     inits_hetero(TLC),
                     params_hetero,
                     path_to_model,
                     n.chains = fitting_options[["nc"]],
                     n.thin = fitting_options[["nt"]],
                     n.iter = fitting_options[["ni"]],
                     n.burnin = fitting_options[["nb"]],
                     parallel=T)


}


#' init hetero used internally only
#'
#'matic rounding
#'@name hinit_hetero
#'@param x TLC value for the focal species
#'@return a list with initial values

inits_hetero <- function(x){
  list(
    tauy=runif(1),
    p=runif(2),
    beta1=runif(1,-0.1,0),
    beta2=runif(1,0.05,0.1),
    TMR=runif(1,0.1,0.18),
    tlc=runif(1,x,40)
  )}
