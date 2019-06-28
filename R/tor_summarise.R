#' tor_summarise
#'
#'The function tor_summarise() provides a comprehensive summary of the output
#'returned by fit_torpor(). Values of mean, 95CI bounds, median and
#'Brooks–Gelman–Rubin criterion (i.e. chain convergence check) of parameters
#'posterior distributions are provided. Additionally, prior posterior overlap
#'values for different parameters are generated. Finally, the mean MR value used
#'to standardize the MR variable is reported.
#'
#'@name tor_summarise
#'@aliases tor_summarise
#'@param mod a fitted model with fit_torpor
#'@return a list of 2. The first element is a tibble with the parameter estimates, the second is the overlap values
#'@export
#'@examples
#'data(test_data)
#'test2 <- tor_fit(MR = test_data[,2],
#'Ta = test_data[, 1],
#'BMR = 29,
#'TLC = 28.8,
#'model = NULL,
#'fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'tor_summarise(mod = test2)

tor_summarise <- function(mod){
out <- list()

params <- c("tauy", "inte", "intc", "intr", "betat", "betac", "Tt", "tlc") ## params of interest

mean <- unlist(mod$mean[params])
CI_97.5 <- unlist(mod$q97.5[params])
median <- unlist(mod$q50[params])
CI_2.5 <- unlist(mod$q2.5[params])
Rhat <- unlist(mod$Rhat[params])

## frame the output in a df
x <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
x$parameter <- rownames(x)
x$parameter[x$parameter == "betat"] <- "betat"
x$parameter[x$parameter == "tlc"] <- "Tm"
x <- x[,c(6,1,2,3,4, 5)]
rownames(x) <- NULL

out$parameter_estimates <- x
############### OVERLAP
out$overlap <- tor_overlap(mod)
return(out)

}
