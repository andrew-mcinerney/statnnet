# Methods for statnnet objects

Methods for statnnet objects

## Usage

``` r
# S3 method for class 'statnnet'
print(x, ...)

# S3 method for class 'statnnet'
coef(object, ...)

# S3 method for class 'statnnet'
vcov(object, ...)

# S3 method for class 'statnnet'
predict(object, newdata, ...)

# S3 method for class 'statnnet'
anova(object, ...)

# S3 method for class 'statnnet'
summary(
  object,
  weights = FALSE,
  level = 0.95,
  effects = c("rowwise", "none", "partial"),
  ...
)
```

## Arguments

- x, object:

  A `statnnet` object.

- ...:

  Additional arguments passed to the underlying method where documented.

- newdata:

  Optional data passed unchanged to the retained `nnet` fit.

- weights:

  Include individual-weight Wald diagnostics in the summary.

- level:

  Confidence level for individual-weight intervals and pointwise effect
  intervals. Curve intervals are not simultaneous confidence bands.

- effects:

  Effect estimand to include. `"rowwise"` (the default) computes the
  average observed-row predictive finite difference by changing the
  focal predictor within each observed covariate row and then averaging.
  `"partial"` computes the average partial-dependence finite-difference
  contrast used by
  [`pce()`](https://andrew-mcinerney.github.io/statnnet/reference/pce.md),
  which additionally averages over the observed focal-predictor values.
  These estimands generally differ when the fitted response contains
  interactions. `"none"` omits effect summaries.

## References

McInerney, A. and Burke, K. (2026). Investigating Statistical Inference
and Covariate Effects in Shallow Neural Networks. arXiv:2311.08139.
[doi:10.48550/arXiv.2311.08139](https://doi.org/10.48550/arXiv.2311.08139)
.

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
summary(model, effects = "none")
#> statnnet summary
#> Formula: mpg ~ weight + origin
#> Response: continuous; architecture: 3-1-1
#> 
#> Effect summaries: omitted
#> 
#> Grouped input-to-hidden Wald summaries
#>    term             model_columns n_weights effective_df wald_chisq
#>  weight                    weight         1    0.9870157  41.431079
#>  origin originEurope, originJapan         2    1.9921566   3.804563
#>       p_value
#>  1.181581e-10
#>  1.482950e-01
anova(model)
#> Grouped input-to-hidden Wald summaries
#>    term             model_columns n_weights effective_df wald_chisq
#>  weight                    weight         1    0.9870157  41.431079
#>  origin originEurope, originJapan         2    1.9921566   3.804563
#>       p_value
#>  1.181581e-10
#>  1.482950e-01
coef(model)
#>       b->h1      i1->h1      i2->h1      i3->h1        b->o       h1->o 
#>  0.61784268  1.15065462 -0.06435239 -0.20420256  2.20067835 -3.68083437 
vcov(model)
#>               b->h1        i1->h1        i2->h1       i3->h1         b->o
#> b->h1   0.030865695 -0.0026179158  0.0014273309  0.002867042  0.041622307
#> i1->h1 -0.002617916  0.0319568323  0.0002843053 -0.001255788 -0.052660319
#> i2->h1  0.001427331  0.0002843053  0.0105621334  0.004943017  0.008104294
#> i3->h1  0.002867042 -0.0012557880  0.0049430173  0.011194215  0.013673241
#> b->o    0.041622307 -0.0526603191  0.0081042937  0.013673241  0.145201363
#> h1->o  -0.031881655  0.0743051691 -0.0081201455 -0.014903946 -0.170531588
#>               h1->o
#> b->h1  -0.031881655
#> i1->h1  0.074305169
#> i2->h1 -0.008120146
#> i3->h1 -0.014903946
#> b->o   -0.170531588
#> h1->o   0.217959163
head(predict(model))
#>         [,1]
#> 1 -0.7108587
#> 2 -0.8551218
#> 3 -0.6533021
#> 4 -0.6506937
#> 5 -0.6645376
#> 6 -1.1917182
```
