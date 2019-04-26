#'toRpoR: Binary mixture model aimed at distinguishing torpid and
#'euthermic metabolic rates
#'
#'This package enables an objective and standardized distinction between torpid
#'and euthermic metabolic rates (MR) measured in steady-state conditions.
#'Furthermore, it provides parameters’ estimations of the relation between
#'ambient temperatures (Ta) and MR in both physiological stages. This package
#'is aimed to support any physiologist working in thermal energetics. More information
#'can be found on the compagnon article Fasel et al. (in prep)
#'
#'This R package is center around the XXX functions which
#'allow you to fit binary mixture model on metabolic rates data using Bayesian inference.
#'
#'Note: This package is not ready for general use yet!
#'
#'@docType package
#'@name toRpoR
#'
#'@section fit_torpor:
#'The function [fit_torpor()] considers the assumed relation between
#'MR and Ta (Speakman & Thomas 2003).In the hypothermic state (torpor) and above some threshold Ta (Tmin),
#'MR follows an exponential curve reflecting the Arrhenius rate enhancing
#'effect of temperature on chemical reactions, whereas below Tmin, it increases
#'linearly with decreasing Ta to maintain a minimal Tb in torpor. In the
#'euthermic state, MR solely increases linearly with decreasing Ta.
#'
#'@section plot_torpor:
#'The function [plot_torpor()] is a wrapper function arround the [fit_torpor()] and [get_prediction()].
#'It uses [fit_torpor()] to fit a binomial mixture model using
#'Bayesian inference and plot the predicted value as well as the raw data.
#'Measures are presented in different colors depending of the metabolic stage
#'and predicted values as well as 95% credibility interval (segmented lines)
#'are presented. This function enable the user to replicates the analysis done in
#'Fasel et al. (in prep).
#'
#'@section get_prediction:
#'The function provides the predicted MR and 95% credible interval boundaries
#'at a defined Ta given a certain model, in normothermic and/or torpid stage.
#'
#'@section get_classification:
#'The function provides the classification of the individual points. It gives
#'the estimated state.
#'
#'@section get_summary:
#'The function gives a summary of the model fit.
NULL
