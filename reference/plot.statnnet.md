# Plot partial covariate effects

Plot partial covariate effects

## Usage

``` r
# S3 method for class 'statnnet'
plot(
  x,
  variable = NULL,
  uncertainty = c("delta", "simulation", "none"),
  level = 0.95,
  by = NULL,
  ...
)
```

## Arguments

- x:

  A `statnnet` object.

- variable:

  One or more predictor names. The default plots all supported
  predictors.

- uncertainty:

  One of `"delta"`, `"simulation"`, or `"none"`.

- level:

  Pointwise confidence level. Curve intervals are not simultaneous
  confidence bands.

- by:

  Optional conditioning predictor. Two PCE curves are drawn at the
  observed mean minus and plus one standard deviation for a continuous
  predictor, or at both fitted levels for a binary predictor.

- ...:

  Additional arguments passed to
  [`pce()`](https://andrew-mcinerney.github.io/statnnet/reference/pce.md).

## Value

The plotted PCE data, invisibly.

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
plot(model, "weight", length_out = 5)
```
