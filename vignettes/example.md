---
title: "An example"
author: "Colin Vullioud"
date: "2021-05-07"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{An example}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---

## How to use torpor ? An example

This vignette presents a suggested workflow of the `torpor` package. This package is firstly aimed at users not familiar with `rjags`, yet interested in investigating the thermoregulatory pattern of their model species, especially when the latter is a heterotherm and when metabolic rate measurements are to be assigned to either torpor or euthermia.

This workflow is based on one of the provided datasets. We go through a complete analysis to illustrate how to use the package functions. At each step, we also present some options for the users and the possible solution to expected problems. Note, however, that this example does not present all aspects of the package.

## Installing the package
As the package is still not available on CRAN, users can download it by calling the following command.

```r
remotes::install_github("vullioud/torpor", build_vignettes = TRUE ,force=TRUE)
```

The program `JAGS` can be downloaded from the website: [https://mcmc-jags.sourceforge.io/](https://mcmc-jags.sourceforge.io/) 


## Fitting a model  
The first step in any analysis is to import the data. In order to fit a model using the function `tor_fit()`, we need a set of measurements of Metabolic rate $M$ and of the ambient temperature $T_a$ at which the measurements have been recorded. These values can be vectors of the same length or two columns of a data frame. In the chosen example on the Tasmanian pygmy possum, *Cercartetus lepidus*, we have 103 metabolic rate values at various $T_a$ by digitalization of a figure (Geiser, 1987).

The data are accessible with the following command: 

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
The lower critical temperature ($T_{lc}$) of the thermoneutral zone ($TNZ$) and the metabolic rate within $TNZ$ ($M_{TNZ}$)  can be estimated by or provided to the `tor_fit()` function. These two values were originally estimated by the author and can be found in the documentation of the dataset (test_data2) for the present example, we will estimate $T_{lc}$ and provide $M_{TNZ}$.
`tor_fit()` represents the core function of the model and should be the first step in any analysis using the torpor package. 

The model is fitted with the following call: 


```r
model <-tor_fit(Ta = test_data2$Ta, 
                M = test_data2$VO2ms,
                Mtnz = 1.8,
                fitting_options = list(ni = 500,
                                            nt = 2,
                                            nb = 200,
                                            nc = 2,
                                            parallel =FALSE), 
                confidence = 0.5)
```

```
## Tlc is being estimated from the data
```

```
## $params
##    parameter   mean CI_2.5 median CI_97.5  Rhat
## 1       tau1  0.148  0.098  0.141   0.235 1.308
## 2       tau2  0.429  0.339  0.422   0.551 1.092
## 3       tau3  0.136  0.092  0.133   0.191 1.007
## 4       inte  7.435  7.168  7.437   7.756 1.055
## 5       intc  0.099  0.075  0.090   0.144 1.995
## 6       intr  0.455 -0.960  0.923   1.203 2.561
## 7      betat -0.199 -0.210 -0.199  -0.190 1.055
## 8      betac  0.065  0.029  0.068   0.093 1.541
## 9         Tt  1.730 -5.362  4.022   5.516 2.559
## 10       TMR  0.114  0.087  0.114   0.146 1.363
## 11        Mr  0.641  0.265  0.638   1.092 1.178
## 12       Tbe 37.355 36.862 37.346  37.798 1.052
## 13       Tbt  2.302 -4.865  4.632   6.155 2.560
## 14       Tlc 28.253 26.837 28.307  29.999 2.689
## 15      Mtnz  1.800     NA  1.800      NA    NA
## 
## $ppo
##   name  ppo
## 1   Mr 34.8
## 2  TMR  4.3
```
Note that for the purpose of this vignette we only run two chains with 500 iterations and 200 burned-in. The results have to be considered carefully (i.e. Checking the adequacy of the model, below). Default settings are more suitable but require more time for the model to run.

The output of the `tor_fit()` unction is a list containing information on the posterior distributions of the investigated parameters as well as on the convergence of the chains and on the prior-posterior distributions overlap of some parameters. These information can also be called by the function  `tor_summarise()`. The rest of the package offers several functions to make more sense of the model output, to check its consistency and to represent it graphically. Researchers familiar with `jagsUI` can develop their analysis from that point.

## Checking the adequacy of the model 

Once the model is fitted, it is recommended to verify the adequacy of the results. This can be achieved via the function `tor_summarise()`. The latter will return a list with the essentials: A dataframe with the mean and median of parameter estimates, the 95% credible interval and the $\hat{R}$ value (i.e. chain convergence estimation). It also reports the parameters’ identifiability.


```r
summary <- tor_summarise(model)
```

The parameter estimates are accessible with: 


```r
summary$params
```

```
##    parameter   mean CI_2.5 median CI_97.5  Rhat
## 1       tau1  0.148  0.098  0.141   0.235 1.308
## 2       tau2  0.429  0.339  0.422   0.551 1.092
## 3       tau3  0.136  0.092  0.133   0.191 1.007
## 4       inte  7.435  7.168  7.437   7.756 1.055
## 5       intc  0.099  0.075  0.090   0.144 1.995
## 6       intr  0.455 -0.960  0.923   1.203 2.561
## 7      betat -0.199 -0.210 -0.199  -0.190 1.055
## 8      betac  0.065  0.029  0.068   0.093 1.541
## 9         Tt  1.730 -5.362  4.022   5.516 2.559
## 10       TMR  0.114  0.087  0.114   0.146 1.363
## 11        Mr  0.641  0.265  0.638   1.092 1.178
## 12       Tbe 37.355 36.862 37.346  37.798 1.052
## 13       Tbt  2.302 -4.865  4.632   6.155 2.560
## 14       Tlc 28.253 26.837 28.307  29.999 2.689
## 15      Mtnz  1.800     NA  1.800      NA    NA
```

### Checking the convergence
If the $\hat{R}$ given in the summary dataframe are larger than 1.1, we recommend to refit the model and increase the number of iterations and burn-ins respectively. This can be done by setting the parameter  `fitting_options = list(ni = , nb = )`. It is also possible to look graphicaly at the convergence using the jagsUI function `traceplot()`. Let’s have a look at the convergence for the parameters betat and betac. The necessary object to run the `jagsUI::traceplot()` function can be accessed through `model$mod_parameter`.


```r
jagsUI::traceplot(model$mod_parameter, c("betat", "betac"))
```

### checking the parameters identifiability 
In addition to a verification of the convergence it is advised to control the identifiability of some parameters. This is done by comparing the prior and the posterior distributions. A warning will be sent to the user for estimated parameters whose PPO are higher than 75% (Fasel et al. in prep.). Parameter identifiability is provided for $T_{lc}$, $M_r$, and $T_{be}$, $TMR$ and $T_{bt}$. The ovelap is also given by the `tor_summarise()` function. Let’s continue with our analysis of *Cercartetus lepidus*:


```r
summary$ppo
```

```
##   name  ppo
## 1   Mr 35.5
## 2  TMR  4.2
```

## Making sense of the model 

Once the basic checks have been done, we can go on with the evaluation of the output. There are two main functions that deal with predictions. `tor_assign()` firstly gives the classification of the raw data based on the predictions of the model. `tor_predict()` further gives the predicted M in torpor and euthermia for a given $T_a$. Finally, we can plot the result using the function `tor_plot()`.

### Getteing the assignation  
TTo look at the assignments of the data the corresponding predicted values we use the function `tor_assign()`,  , which will assign the measured metabolic values to either torpor or euthermia  and returns a dataframe with the measured $M$, the measured $T_a$, the predicted $M$ and the assignment. (predicted_state).


```r
assignment <- tor_assign(tor_obj = model)
 
head(assignment)
```

```
##   measured_M measured_Ta assignment predicted_M
## 1  0.3643880     2.81631     Torpor   0.3571705
## 2  0.2222220     3.32812     Torpor   0.2563818
## 3  0.5080130     3.37791     Torpor   0.2466737
## 4  0.4407550     3.94446     Torpor   0.1592819
## 5  7.0744500     4.55155  Euthermia   6.5304441
## 6  0.0538837     4.60326     Torpor   0.1324550
```

### Prediction 
The `tor_predict()`  function is slightly different as it takes a vector of $T_a$ as input and returns the predicted  $M$ both in torpor and in euthermia and the 95% credible interval. For example,  let’s see what are the predicted $M$ at $T_a$ 22°!


```r
prediction <- tor_predict(tor_obj = model, Ta = 22)
head(prediction)
```

```
##   Ta assignment      pred    upr_95    lwr_95
## 1 22     Torpor 0.4103882 0.6136446 0.2177884
## 2 22  Euthermia 3.0558678 3.1269018 2.9960147
```

### plotting the data and predicted values 
Finally a built-in function allows plot the results. The user can modulate labels `xlab` and `ylab` and the colors with `col_torp`, `col_eut` and `col_Mtnz` and can save the plot using the arguments `pdf = TRUE`. The plot_type argument also allows to choose between `base` (default) and `ggplot`.


```r
 tor_plot(tor_obj = model, ylab = "MR")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)





