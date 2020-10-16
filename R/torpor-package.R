#'torpor: Mixture model aimed at assigning metabolic rate measurements (M) to torpor and euthermia in heterothermic endotherms.
#'#'
#'#'This package estimate the assignment of M to torpor or euthermia.
#'It uses the variation in M measured during euthermic rest and torpor at
#'different ambient temperatures (Ta) to estimate the lower critical
#'temperature (Tlc) of the thermoneutral zone (TNZ) and determine physiological
#'state membership using mixture models. In addition, this package allows
#'the further prediction of M during rest and torpor along Ta, including resting
#'metabolic rate within the TNZ.
#'
#'This package is aimed to support any physiologist working in thermal energetics.
#'More information can be found on the companion article Fasel et al. (in prep)
#'and in the vignettes.
#'
#'This package is center around the [tor_fit()] function which allows to fit
#'mixture models on metabolic rates data using Bayesian inference.
#'@docType package
#'@name torpor
#'
#'@section tor_fit:
#'The function [tor_fit()] considers the relation between metabolic rate (M) and
#'ambient temperature (Ta) assumed by the Scholander-Irving model and its later extensions.
#'
#'Resting M measured within the thermoneutral zone (TNZ) is independent of Ta.
#'This rate is hereafter referred to as MTNZ, although it would correspond to
#'the basal metabolic rate (BMR) provided that the specific criteria for the
#'BMR are met (see Fasel et al. in prep.). Below the lower critical temperature
#'of TNZ (Tlc), M of euthermic animals increases linearly with decreasing Ta.
#'M of torpid animals increases linearly with decreasing Ta to maintain a minimal
#'body temperature below some threshold ambient temperature (Tt). This state is
#'usually referred to as "regulated torpor". Between Tt and Tlc, M of torpid
#'animals follows an exponential curve. In this Ta range, torpor is referred to
#'as "conforming torpor".
#'
#'@section tor_plot:
#'The function [tor_plot()] is a wrapper function around the [tor_fit())] and [tor_predict()].
#'
#'It uses [tor_fit()] to fit a mixture model using
#'Bayesian inference and plot the predicted value as well as the raw data.
#'Measures are presented in different colors depending of the metabolic state.
#'Predicted values as well as 95% credibility interval (segmented lines)
#'are also presented. This function enables the user to replicates the analysis done in
#'Fasel et al. (in prep).
#'
#'@section tor_predict:
#'The function provides the predicted M and 95% credible interval boundaries
#'at a defined Ta given a certain model, in euthermic and/or torpid state.
#'
#'@section tor_classify:
#'The function classify the individual points according to their estimated state.
#'
#'@section tor_summarise:
#'The function gives a summary statistics of the model fit.
#'NULL
