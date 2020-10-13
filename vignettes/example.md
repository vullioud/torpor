---
title: "An example"
author: "Colin Vullioud"
date: "2020-10-13"
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
## Warning in estimate_tlc_mtnz(Ta = Ta, M = Y, fitting_options = fitting_options):
## Mtnz computed on less than 10 points
```

```
## Warning in eval(family$initialize): non-integer #successes in a binomial glm!
```

```
## Warning in get_parameters(tor_obj): Values for Mtnz are directly computed from
## data points
```

```
## $params
##    parameter   mean CI_2.5 median CI_97.5  Rhat
## 1       tau1  0.101  0.076  0.098   0.140 1.005
## 2       tau2  0.395  0.318  0.391   0.505 1.006
## 3       inte  7.391  7.146  7.398   7.680 0.998
## 4       intc  0.033  0.014  0.031   0.061 1.079
## 5       intr  0.945  0.811  0.922   1.194 1.011
## 6      betat -0.196 -0.206 -0.196  -0.187 0.998
## 7      betac  0.117  0.084  0.116   0.148 1.043
## 8         Tt  4.552  3.961  4.457   5.725 1.025
## 9        TMR  0.055  0.028  0.054   0.091 1.039
## 10        MR  0.894  0.688  0.888   1.121 1.020
## 11       Tbe 37.769 37.326 37.752  38.171 0.998
## 12       Tbt  4.832  4.171  4.717   6.066 1.012
## 13       Tlc 28.937 26.743 28.799  31.682 1.029
## 14      MTNZ  1.755  1.567  1.746   1.942    NA
## 
## $ppo
##   name  ppo
## 1  Tlc 35.5
## 2   MR 26.0
## 3  Tbe  5.9
## 4  TMR  5.7
## 5  Tbt  6.5
```
Note that for the purpose of this example we only run 5000 iterations with 3000 
burned-in. The results are thus not trusthworthy.

The output of the `tor_fit()` function is a list of class `tor_obj`. The rest of the package offers several functions to manipulate the model output, to check its consistency and represent it graphically.

## Checking the robustness of the model 

Once the model is fitted, it is recommended to have an overlook at the results. This can be achieved via the function `tor_summarise()`. The latter will return a list with the essentials: A data.frame with the mean and median of parameter estimates, the 95% credible interval and the $\hat{R}$ value (i.e. chain convergence estimation). It also reports the parameters’ identifiability. 


```r
summary <- tor_summarise(model)
```

```
## Warning in get_parameters(tor_obj): Values for Mtnz are directly computed from
## data points
```

The parameter estimates are accessible with: 


```r
summary$params
```

```
##    parameter   mean CI_2.5 median CI_97.5  Rhat
## 1       tau1  0.101  0.076  0.098   0.140 1.005
## 2       tau2  0.395  0.318  0.391   0.505 1.006
## 3       inte  7.391  7.146  7.398   7.680 0.998
## 4       intc  0.033  0.014  0.031   0.061 1.079
## 5       intr  0.945  0.811  0.922   1.194 1.011
## 6      betat -0.196 -0.206 -0.196  -0.187 0.998
## 7      betac  0.117  0.084  0.116   0.148 1.043
## 8         Tt  4.552  3.961  4.457   5.725 1.025
## 9        TMR  0.055  0.028  0.054   0.091 1.039
## 10        MR  0.894  0.688  0.888   1.121 1.020
## 11       Tbe 37.769 37.326 37.752  38.171 0.998
## 12       Tbt  4.832  4.171  4.717   6.066 1.012
## 13       Tlc 28.937 26.743 28.799  31.682 1.029
## 14      MTNZ  1.755  1.567  1.746   1.942    NA
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
## 1  Tlc 39.3
## 2   MR 25.4
## 3  Tbe  6.6
## 4  TMR  5.5
## 5  Tbt  6.6
```

### Making sens of the model 

Once the basic checks have been done, we can go on with the evaluation of the output. There are two main functions that deal with predictions. `tor_classify()` firstly gives the classification of the raw data based on the predictions of the model. `tor_predict()` further gives the predicted $M$ in torpor and euthermy for a given $T_a$. Finally, we can plot the result using the function `plot_torpor()`.


### Getteing the classification 

To look at the classification of the data the corresponding predicted values we use the function `tor_classify()`, which will assign the measured metabolic values to either torpor or euthermy and returns a dataframe with the measured $M$, the measured $T_a$, the predicted $M$ and the classification (predicted_state).


```r
classification <- tor_classify(tor_obj = model)
 
head(classification)
```

```
##   measured_M measured_Ta classification predicted_M
## 1  0.3643880     2.81631         Torpor  0.37152688
## 2  0.2222220     3.32812         Torpor  0.27186752
## 3  0.5080130     3.37791         Torpor  0.26234042
## 4  0.4407550     3.94446         Torpor  0.15135214
## 5  7.0744500     4.55155       Euthermy  6.50614920
## 6  0.0538837     4.60326         Torpor  0.06241977
```



### Prediction 
The `tor_predict()` function is slightly different as it takes a vector of $T_a$ as input and return the predicted $M$ both in torpor and in euthermy and the 95% credible interval. For example, let's see what are the predicted $M$ at Ta 22°C. 


```r
prediction <- tor_predict(tor_obj = model, Ta = 22)
head(prediction)
```

```
##   Ta    group     pred    upr_95    lwr_95
## 1 22   Torpor 0.404545 0.4742284 0.3305445
## 2 22 Euthermy 3.086869 3.1534071 3.0273552
```

### plotting the data and predicted values 

Finally a built-in function allows plot the results. The user can specifie some parameters and can save the plot using the arguments `pdf = TRUE`. The `plot_type` argument also allows to choose between `ggplot` and `base-R`


```r
plot <- tor_plot(tor_obj = model, ylab = "MR")
```

```
## Warning in tor_predict(tor_obj, seq(Tlimlo, Tlimup, length = 1000)): Tuc is not
## considered: Mtnz is calculated independently of Ta above Tlc
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

```r
plot 
```

```
## NULL
```





