# Select a parsimonious nnet model

`select_nnet()` implements the information-criterion procedure of
McInerney and Burke. It first selects the hidden-layer size, applies
backward input selection, and then alternates one-step changes to the
hidden and input layers. Each candidate is fitted from several random
initialisations using
[`nnet::nnet()`](https://rdrr.io/pkg/nnet/man/nnet.html).

## Usage

``` r
select_nnet(x, ...)

# S3 method for class 'formula'
select_nnet(
  x,
  data,
  ...,
  max_size,
  n_starts = 5L,
  criterion = c("BIC", "AIC", "AICc"),
  start_range = c(-3, 3),
  maxit = 1000L,
  seed = NULL,
  progress = TRUE
)

# Default S3 method
select_nnet(
  x,
  y,
  ...,
  max_size,
  n_starts = 5L,
  criterion = c("BIC", "AIC", "AICc"),
  start_range = c(-3, 3),
  maxit = 1000L,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- x:

  A model formula or a numeric design matrix.

- ...:

  Arguments passed to methods and then to
  [`nnet::nnet()`](https://rdrr.io/pkg/nnet/man/nnet.html). Arguments
  that would change the supported architecture, response model, initial
  weights, or decay value are not allowed.

- data:

  A data frame used with the formula interface.

- max_size:

  Maximum number of hidden units considered in the initial hidden-layer
  search.

- n_starts:

  Number of random initialisations fitted for each candidate.

- criterion:

  Information criterion used to compare candidates.

- start_range:

  Numeric vector of length two giving the lower and upper bounds for
  uniformly generated starting weights.

- maxit:

  Maximum number of `nnet` optimiser iterations per start.

- seed:

  Optional random seed for reproducible starting weights.

- progress:

  Whether to print progress as candidate models are fitted.

- y:

  A numeric response used with the matrix interface.

## Value

An object of class `"statnnet_selection"` containing the selected native
`nnet` fit, selected terms or columns, the full search path, and
start-level convergence diagnostics.

## Details

Selection is available for unpenalised Gaussian regression. With a
formula, complete right-hand-side terms are selected together: all
columns belonging to a factor or transformed term enter or leave as a
group. With a matrix, columns are selected individually.

## References

McInerney, A. and Burke, K. (2025). A statistical modelling approach to
feedforward neural network model selection. *Statistical Modelling*,
25(4), 323–342.
[doi:10.1177/1471082X241258261](https://doi.org/10.1177/1471082X241258261)
.

## Examples

``` r
set.seed(1)
example_data <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
example_data$y <- 1 + 2 * example_data$x1 + rnorm(60, sd = 0.4)

selection <- select_nnet(
  y ~ x1 + x2,
  data = example_data,
  max_size = 1,
  n_starts = 2,
  maxit = 500,
  seed = 2,
  progress = FALSE
)
selection
#> statnnet model selection
#> Call:
#> select_nnet(x = y ~ x1 + x2, data = example_data, max_size = 1, 
#>     n_starts = 2, maxit = 500, seed = 2, progress = FALSE)
#> Selected model: 1 input column, 1 hidden unit, 4 network weights
#> BIC: 81.590; 5 candidate models; 10 converged starts
#> Selected terms: x1
#> Excluded terms: x2
```
