---
title: "An example"
author: "Colin Vullioud"
date: "2020-12-19"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{An example}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---

## How to use torpor ? An example

This vignette presents a suggested workflow of the `torpor` package. This package is firstly aimed at users not familiar with `rjags`, yet interested in investigating the thermoregulatory pattern of their model species, especially when the latter is a heterotherm and when metabolic rate measurements are to be assigned to either torpor or euthermia.

This workflow is based on one of the provided datasets. We go through a complete analysis to illustrate how to use the package functions. At each step we also present the options for the users and the possible solution to expected problems. Note, however, that this example does not present all aspects of the package. 

### Fitting a model  
The first step in any analysis is to get data. In order to fit a model using the function `tor_fit()` we need a set of measurements of Metabolic rate $M$ and the ambient temperature $T_a$ at which the measurements have been made. These values can be vectors of the same length or two columns of a data.frame. In the chosen example on the Tasmanian pygmy possum, *Cercartetus lepidus*, we have 103 metabolic rate values at various $T_a$ by digitalization of a figure (Geiser, 1987).

The data are accessible with the following command. 


```r
library(torpor)
data(test_data2)
str(test_data2)
```

```
## 'data.frame':	103 obs. of  2 variables:
##  $ Ta   : num  2.82 3.01 3.09 3.31 3.33 ...
##  $ VO2ms: num  0.364 NA NA NA 0.222 ...
```

The lower critical temperature (i.e., lower limit of the thermoneutral zone, $T_{lc}$ and the basal metabolic rate $Mtnz$ are also needed in order to fit the model with the `tor_fit()` function. These two values can be found in the documentation of the dataset (?test_data2) for the present example and should be provided by the researcher when using her/his own data. 
`tor_fit()` represents the core function of the model and should be the first step in any analysis using the torpor package. Researchers not familiar with `R` and only interested in a visual representation of their data can skip this step and use the `tor_plot()` function directly. This option is however not recommended since problems might remain unnoticed. 

The model is fitted with the following call: 


```r
model <- tor_fit(M = test_data2[, "VO2ms"], Ta = test_data2[, "Ta"],
                    fitting_options = list(ni = 5000,
                                           nt = 10,
                                           nb = 3000,
                                           nc = 2))
```

```
## Mtnz and Tlc are being estimated from the data
```

```
## Warning in estimate_Tlc_Mtnz(Ta = Ta, M = Y, Mtnz = Mtnz, fitting_options =
## fitting_options): Mtnz computed on less than 10 points
```

```
## Warning in eval(family$initialize): non-integer #successes in a binomial glm!
```

```
## $params
##    parameter   mean CI_2.5 median CI_97.5  Rhat
## 1       tau1  0.129  0.092  0.126   0.199 1.033
## 2       tau2  0.444  0.350  0.439   0.562 1.002
## 3       inte  7.261  6.981  7.259   7.558 1.007
## 4       intc  0.081  0.068  0.078   0.112 1.020
## 5       intr  0.895  0.076  0.921   1.193 1.038
## 6      betat -0.186 -0.196 -0.186  -0.176 1.007
## 7      betac  0.075  0.055  0.077   0.090 1.007
## 8         Tt  4.212 -0.113  4.347   5.870 1.037
## 9        TMR  0.112  0.099  0.109   0.139 1.011
## 10        MR  0.759  0.508  0.762   1.015 1.009
## 11       Tbe 39.092 38.602 39.086  39.587 1.007
## 12       Tbt  4.814  0.416  4.946   6.497 1.036
## 13       Tlc 29.667 28.053 29.662  31.326 1.013
## 14      Mtnz  1.750  1.528  1.701   1.972    NA
## 
## $ppo
##   name  ppo
## 1  Tlc 25.6
## 2   MR 28.3
## 3  Tbe  7.4
## 4  TMR  4.5
## 5  Tbt 11.7
```
Note that for the purpose of this example we only run 5000 iterations with 3000 
burned-in. The results are thus not trusthworthy.

The output of the `tor_fit()` function is a list of class `tor_obj`. The rest of the package offers several functions to manipulate the model output, to check its consistency and represent it graphically.

## Checking the robustness of the model 

Once the model is fitted, it is recommended to have an overlook at the results. This can be achieved via the function `tor_summarise()`. The latter will return a list with the essentials: A data.frame with the mean and median of parameter estimates, the 95% credible interval and the $\hat{R}$ value (i.e. chain convergence estimation). It also reports the parameters’ identifiability. 


```r
summary <- tor_summarise(model)
```

The parameter estimates are accessible with: 


```r
summary$params
```

```
##    parameter   mean CI_2.5 median CI_97.5  Rhat
## 1       tau1  0.129  0.092  0.126   0.199 1.033
## 2       tau2  0.444  0.350  0.439   0.562 1.002
## 3       inte  7.261  6.981  7.259   7.558 1.007
## 4       intc  0.081  0.068  0.078   0.112 1.020
## 5       intr  0.895  0.076  0.921   1.193 1.038
## 6      betat -0.186 -0.196 -0.186  -0.176 1.007
## 7      betac  0.075  0.055  0.077   0.090 1.007
## 8         Tt  4.212 -0.113  4.347   5.870 1.037
## 9        TMR  0.112  0.099  0.109   0.139 1.011
## 10        MR  0.759  0.508  0.762   1.015 1.009
## 11       Tbe 39.092 38.602 39.086  39.587 1.007
## 12       Tbt  4.814  0.416  4.946   6.497 1.036
## 13       Tlc 29.667 28.053 29.662  31.326 1.013
## 14      Mtnz  1.750  1.528  1.701   1.972    NA
```

### Checking the convergence
If the $\hat{R}$ values, given in the summary data.frame are larger than 1.1, we recommend to refit the model and increase the number of iterations and burn-ins respectively. This can be done by setting the parameter `fitting_options = list(ni = , nb = )`

### checking the parameters identifiability 
In addition to check of the convergence it is advised to control the identifiability of the parameters. This is done by comparing the prior and the posterior distributions. If the overlap is greater than 0.3 the parameter should be considered as not identifiable and a warning will be sent to the user. 

The ovelap is also given by the `tor_summarise()` function. Let's continue with our analysis of *Cercartetus lepidus*: 


```r
summary$ppo
```

```
##   name  ppo
## 1  Tlc 26.9
## 2   MR 25.2
## 3  Tbe  6.7
## 4  TMR  3.6
## 5  Tbt 11.2
```

### Making sens of the model 

Once the basic checks have been done, we can go on with the evaluation of the output. There are two main functions that deal with predictions. `tor_classify()` firstly gives the classification of the raw data based on the predictions of the model. `tor_predict()` further gives the predicted $M$ in torpor and euthermy for a given $T_a$. Finally, we can plot the result using the function `plot_torpor()`.


### Getteing the classification 

To look at the classification of the data the corresponding predicted values we use the function `tor_classify()`, which will assign the measured metabolic values to either torpor or euthermy and returns a dataframe with the measured $M$, the measured $T_a$, the predicted $M$ and the classification (predicted_state).







