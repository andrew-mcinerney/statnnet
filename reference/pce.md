# Partial covariate effects for a fitted statnnet model

Computes finite differences of the partial-dependence function. The
results describe the fitted response surface and should not be
interpreted as causal effects. For correlated predictors, some evaluated
combinations may have little support in the data.

## Usage

``` r
pce(
  object,
  variable,
  values = NULL,
  d = NULL,
  length_out = 101L,
  type = c("curve", "average"),
  uncertainty = c("delta", "simulation", "none"),
  level = 0.95,
  nsim = 1000L,
  seed = NULL,
  by = NULL
)
```

## Arguments

- object:

  A `statnnet` object.

- variable:

  Name of an original predictor.

- values:

  Values at which to evaluate a continuous predictor. By default, an
  evenly spaced grid keeps the baseline and finite-difference comparison
  within the observed range. User-supplied values are used unchanged and
  may extrapolate.

- d:

  Finite-difference step for a continuous predictor. The default is one
  observed standard deviation. Binary predictors are always compared
  from 0 to 1 and factors against their first level.

- length_out:

  Number of grid points for a continuous predictor.

- type:

  `"curve"` for evaluations over `values`, or `"average"` for a point
  summary averaged over observed predictor values. Average effects may
  compare values beyond the observed range near a boundary.

- uncertainty:

  One of `"delta"`, `"simulation"`, or `"none"`.

- level:

  Pointwise confidence level. Curve intervals are not simultaneous
  confidence bands.

- nsim:

  Number of simulation draws when `uncertainty = "simulation"`.

- seed:

  Optional simulation seed.

- by:

  Optional name of a second predictor at whose values the PCE is
  evaluated. Continuous predictors use their observed mean minus and
  plus one standard deviation. Binary predictors use both fitted levels.

## Value

A data frame of class `"statnnet_pce"`.

## Examples

``` r
data("auto_mpg", package = "statnnet")
example_data <- auto_mpg
example_data$mpg <- as.numeric(scale(example_data$mpg))
example_data$weight <- as.numeric(scale(example_data$weight))
set.seed(3)
fit <- nnet::nnet(
  mpg ~ weight + origin,
  data = example_data,
  size = 1,
  linout = TRUE,
  decay = 0.1,
  Hess = TRUE,
  maxit = 500,
  trace = FALSE
)
model <- statnnet(fit, data = example_data, response = "continuous")
pce(model, "weight", length_out = 5)
#>   variable contrast       value step   estimate  std_error   conf_low
#> 1   weight     <NA> -1.60652233    1 -0.9190229 0.10378270 -1.1224333
#> 2   weight     <NA> -0.81843976    1 -1.0197353 0.06368578 -1.1445572
#> 3   weight     <NA> -0.03035718    1 -0.7852639 0.04444044 -0.8723656
#> 4   weight     <NA>  0.75772540    1 -0.4496493 0.06474863 -0.5765543
#> 5   weight     <NA>  1.54580798    1 -0.2142032 0.06056191 -0.3329024
#>     conf_high
#> 1 -0.71561257
#> 2 -0.89491351
#> 3 -0.69816226
#> 4 -0.32274432
#> 5 -0.09550409
```
