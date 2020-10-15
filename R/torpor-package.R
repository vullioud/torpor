#'torpor: Binary mixture model aimed at distinguishing torpid and
#'euthermic metabolic rates
#'
#'This package enables an objective and standardized distinction between torpid
#'and euthermic metabolic rates (M) measured in steady-state conditions.
#'Furthermore, it provides parameters’ estimations of the relation between
#'ambient temperatures (Ta) and MR in both physiological stages. This package
#'is aimed to support any physiologist working in thermal energetics. More information
#'can be found on the compagnon article Fasel et al. (in prep) and in the vignettes.
#'
#'This package is center around the [tor_fit()] function which
#'allows to fit binary mixture model on metabolic rates data using Bayesian inference.
#'
#'@docType package
#'@name torpor
#'
#'@section tor_fit:
#'The function [tor_fit()] considers the assumed relation between
#'metabolic rate (M) and ambient temperature (Ta) (Speakman & Thomas 2003).In the hypothermic state (torpor) and above some threshold Ta (Tmin),
#'M follows an exponential curve reflecting the Arrhenius rate enhancing
#'effect of temperature on chemical reactions, whereas below Tmin, it increases
#'linearly with decreasing Ta to maintain a minimal Tb in torpor. In the
#'euthermic state, M solely increases linearly with decreasing Ta.
#'
#'@section tor_plot:
#'The function [tor_plot()] is a wrapper function around the [tor_fit())] and [tor_predict()].
#'It uses [tor_fit()] to fit a binomial mixture model using
#'Bayesian inference and plot the predicted value as well as the raw data.
#'Measures are presented in different colors depending of the metabolic state.
#'Predicted values as well as 95% credibility interval (segmented lines)
#'are also presented. This function enable the user to replicates the analysis done in
#'Fasel et al. (in prep).
#'
#'@section tor_predict:
#'The function provides the predicted M and 95% credible interval boundaries
#'at a defined Ta given a certain model, in normothermic and/or torpid stage.
#'
#'@section tor_classify:
#'The function classify the individual points according to their estimated state.
#'
#'@section tor_summarise:
#'The function gives a summary statistics of the model fit.
NULL
