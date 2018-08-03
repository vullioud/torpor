#' toRpoR: a package to fit mixture models on bats
#'
#' This R package is center around 3 main functions:
#' the fit_torpor() functions allow you to fit gives the tools to reproduce the main statistical analysis of
#' the article: Super Super
#'
#' This package is not conceived for general use yet!
#'
#' In the examples below, we provide the workflow leading the main results
#' presented in the paper.
#'
#'
#' NOTE: IT is a mess
#' @docType package
#' @name toRpoR
#'
#' @section fit_torpor():
#' The fit_torpor() function fit a mixture model for 2 pop using jags.
#'
#' @section fit_and_plot():
#' The fit_and_plot() isa wrapper function around the fit_torpor() function. It
#' reproduce the plot presented in the article.
#' @section get_prediction():
#' The get_prediction() functions allow you to access the prediction of the
#' fitted model and the 95 credibility intervals.
#' @section get_Q10():
#' The get_Q10() function gives you the Q10 values for a given dataset.
#'
NULL
