
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" height="139" />

# statnnet

<!-- badges: start -->

[![R-CMD-check](https://github.com/andrew-mcinerney/statnnet/workflows/R-CMD-check/badge.svg)](https://github.com/andrew-mcinerney/statnnet/actions)
<!-- badges: end -->

## Installation

You can install the development version of statnnet from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("andrew-mcinerney/statnnet")
```

## statnnet()

The primary function in this package is `statnnet()`. It creates a
statistically-based version of an existing `nnet` object.

``` r
library(statnnet)
stnn <- selectnn(nnet, X)
```

A useful summary table can be generated using

``` r
summary(stnn)
```

and covariate-effect plots can be created using

``` r
plot(stnn, conf_int = TRUE)
```

More information about these functions and their arguments can be found
in the function documentation.
