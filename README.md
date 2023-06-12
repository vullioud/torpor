
<!-- README.md is generated from README.Rmd. Please edit that file -->

# torpor

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Build
Status](https://travis-ci.org/vullioud/torpor.svg?branch=master)](https://travis-ci.org/vullioud/torpor)
[![Codecov test
coverage](https://codecov.io/gh/vullioud/torpor/branch/master/graph/badge.svg)](https://codecov.io/gh/vullioud/torpor?branch=master)

## Welcome to torpor

This is the repository of torpor an `R` (<https://www.r-project.org/>)
package aiming at an objective and standardized distinction between
torpid and euthermic metabolic rates (M) measured in steady-state
conditions. Note: This package is not aimed for use without careful
researcher’s attention.

A stable version of the package will be soon available on CRAN. However,
to allow greater flexibility and interactivity at this early stage of
the project, the package is only available on github at the moment. Any
bug report or request features would therefore be appreciated.

## Installation

Torpor is using `JAGS` (Just Another Gibbs Sampler) in the background
and thus you should make sure to have
[JAGS](http://mcmc-jags.sourceforge.net) installed on your machine. Once
jags is properly installed and working you can installed torpor in R by
typing:

``` r
library(remotes)
remotes::install_github("vullioud/torpor", build_vignettes = TRUE)
```

## Where can you find more about torpor?

You can learn more about the goal and the theory behind the package in
the companion article (Fasel et al., Biol Open 15 April 2022; 11 (4):
bio059064. doi: <https://doi.org/10.1242/bio.059064>). And you will find
an example in the vignettes of the package. You can access it directly
from R by typing:

``` r
browseVignettes("torpor")
```

## How can you help ?

Torpor is still in its infancy and despite numerous testing it is
possible that bugs appear. You can help us by reporting bugs or request
new features.
