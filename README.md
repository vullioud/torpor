
<!-- README.md is generated from README.Rmd. Please edit that file -->

# torpor

<!-- badges: start -->

<!-- badges: end -->

## Welcome to torpor

This is the repository of torpor an `R` (<https://www.r-project.org/>)
package aiming at an objective and standardized distinction between
torpid and euthermic metabolic rates (MR) measured in steady-state
conditions. Note: This package is not aimed for general use without
careful researcher’s attention.

## Installation

ToRpoR is using jags in the background and thus you should make sure to
have [jags](http://mcmc-jags.sourceforge.net) installed on your machine.
Once jags is properly installed and working you can installed torpor in
R by typing:

``` r
library(remotes)
remotes::install_github("vullioud/torpor", build_opts =  c("--build_vignettes"),force=TRUE)
#> Downloading GitHub repo vullioud/torpor@master
#> 
#>   
   Warning: unknown option ‘--build_vignettes’
#> 
  
   checking for file ‘/tmp/RtmpFRRx4y/remotes6a303b7f64d8/vullioud-torpor-251f7f8/DESCRIPTION’ ...
  
   checking for file ‘/tmp/RtmpFRRx4y/remotes6a303b7f64d8/vullioud-torpor-251f7f8/DESCRIPTION’ ... 
  
✔  checking for file ‘/tmp/RtmpFRRx4y/remotes6a303b7f64d8/vullioud-torpor-251f7f8/DESCRIPTION’
#> 
  
─  preparing ‘torpor’:
#> 
  
✔  checking DESCRIPTION meta-information
#> 
  
─  installing the package to build vignettes
#> 
  
   creating vignettes ...
  
✔  creating vignettes (40.4s)
#> 
  
─  checking for LF line-endings in source and make files and shell scripts (452ms)
#> 
  
─  checking for empty or unneeded directories
#> ─  looking to see if a ‘data/datalist’ file should be added
#> 
  
─  building ‘torpor_0.1.1.tar.gz’
#> 
  
   
#> 
#> Installing package into '/home/colin/R/x86_64-pc-linux-gnu-library/3.6'
#> (as 'lib' is unspecified)
```

The package will be submitted to cran in the forseeable future.

## Where can I find more about toRpoR?

You can learn more about the goal and the theory behind the package in
the companion article (Fassel et al., in prep).
