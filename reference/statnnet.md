# Augment a fitted nnet model with statistical summaries

`statnnet()` adds covariance estimates, Wald summaries, and
covariate-effect calculations to a supported model fitted by
[`nnet::nnet()`](https://rdrr.io/pkg/nnet/man/nnet.html). The original
fit is retained. A model selected by
[`select_nnet()`](https://andrew-mcinerney.github.io/statnnet/reference/select_nnet.md)
can also be supplied.

## Usage

``` r
statnnet(object, ...)

# Default S3 method
statnnet(object, ...)

# S3 method for class 'statnnet_selection'
statnnet(object, covariance_tol = 1e-10, objective_tol = 1e-06, ...)

# S3 method for class 'nnet'
statnnet(
  object,
  formula = NULL,
  data = NULL,
  response = NULL,
  covariance_tol = 1e-10,
  objective_tol = 1e-06,
  ...
)
```

## Arguments

- object:

  A fitted object inheriting from `"nnet"`, or a `"statnnet_selection"`
  object returned by
  [`select_nnet()`](https://andrew-mcinerney.github.io/statnnet/reference/select_nnet.md).

- ...:

  Reserved for methods.

- covariance_tol:

  Relative numerical tolerance used to decide whether the penalised
  information matrix is reliably invertible.

- objective_tol:

  Relative tolerance for validating the reconstructed fitting criterion
  against `object$value`.

- formula:

  The formula used to fit `object`. This may be omitted for an
  `nnet.formula` object.

- data:

  The original training data.

- response:

  Either `"continuous"` or `"binary"`. When `NULL`, the response type is
  inferred from the fitted `nnet` object.

## Value

An object of class `"statnnet"` containing the original fit,
reconstructed model information, covariance calculations, and
diagnostics.

## Details

Covariance calculations are returned only when the fitted Hessian passes
the numerical checks. Wald tests and delta-method intervals are still
approximations conditional on the fitted architecture and local
solution. A warning is given when the hidden layer is larger than the
fitted input layer. Binary-response summaries should be treated as
exploratory because their frequentist calibration has not been studied
by the package authors.

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
model
#> statnnet: statistical summary of a fitted nnet model
#> Formula: mpg ~ weight + origin
#> Architecture: 3-1-1; 6 weights; n = 392
#> Response: continuous; decay: 0.1
#> nnet convergence code: 0
#> Covariance: available (reciprocal condition number 0.000415)
```
