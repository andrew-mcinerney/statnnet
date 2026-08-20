
<!-- README.md is generated from README.Rmd. Please edit that file -->

# statnnet

`statnnet` provides regression-style statistical summaries for small,
single-hidden-layer neural networks fitted with `nnet`. It retains the
original fit: it is not a neural-network optimiser or model-selection
package.

## Installation

``` r
# install.packages("pak")
pak::pak("andrew-mcinerney/statnnet")
```

## Basic workflow

``` r
library(nnet)
library(statnnet)

set.seed(2026)
fit <- nnet(
  medv ~ rm + lstat + chas,
  data = Boston,
  size = 1,
  linout = TRUE,
  decay = 0.1,
  Hess = TRUE,
  maxit = 1000,
  trace = FALSE
)

model <- statnnet(
  fit,
  formula = medv ~ rm + lstat + chas,
  data = Boston,
  response = "continuous"
)

model
#> statnnet: statistical summary of a fitted nnet model
#> Formula: medv ~ rm + lstat + chas
#> Architecture: 3-1-1; 6 weights; n = 506
#> Response: continuous; decay: 0.1
#> nnet convergence code: 0
#> Covariance: available (reciprocal condition number 0.000367)
anova(model)
#> Grouped input-to-hidden Wald summaries
#>   term model_columns n_weights effective_df wald_chisq      p_value status
#>     rm            rm         1    0.9982582   62.67532 2.425982e-15     ok
#>  lstat         lstat         1    0.9965359   99.36074 2.083162e-23     ok
#>   chas          chas         1    0.9995955   12.92108 3.246427e-04     ok
```

`summary(model)` returns structured average partial covariate effects
and grouped Wald-type summaries. `pce()` and `plot()` provide
finite-difference partial covariate effects with delta-method or
simulation uncertainty.

``` r
summary(model)
pce(model, "rm")
plot(model, "rm")
plotnn(model, annotate = "p_value")
```

These summaries describe the fitted response function. They are
conditional on the chosen network and are not post-selection or causal
inference.
