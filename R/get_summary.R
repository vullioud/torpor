#' get_summary
#'
#' This function return a comprehensive summary of the model fit. It returns the
#' estimated parameters of the model. The check of overlap
#'@name get_summary
#'@aliases get_summary
#'@param mod a fitted model with fit_torpor
#'@return a list
#'@import purrr
#'@import overlapping
#'@import dplyr
#'@export
#'@examples
#'data(test_data)
#'test2 <- fit_torpor(MR = test_data[,2],
#'Ta = test_data[, 1],
#'BMR = 29,
#'TLC = 28.8,
#'Model = NULL,
#'fitting_options = list(nc = 1))
#'get_summary(mod = test2)

get_summary <- function(mod){
out <- list()

params <- c("tauy", "inte", "intc", "intr", "betat", "betac", "Tt", "tlc")

mean <- unlist(mod$mean[params])
CI_97.5 <- unlist(mod$q97.5[params])
median <- unlist(mod$q50[params])
CI_2.5 <- unlist(mod$q2.5[params])
Rhat <- unlist(mod$Rhat[params])

## frame the output in a df
x <- as.data.frame(cbind(mean, CI_2.5, median, CI_97.5, Rhat))
x$parameter <- rownames(x)
x <- x[,c(6,1,2,3,4, 5)]
rownames(x) <- NULL

out$parameter_estimates <- x
############### OVERLAP
out$overlap <- check_overlap(mod)
out$Ym <- mod$mean$Ym
##############
return(out)

}
