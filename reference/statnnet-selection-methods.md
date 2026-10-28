# Methods for statnnet model-selection results

Methods for statnnet model-selection results

## Usage

``` r
# S3 method for class 'statnnet_selection'
print(x, ...)

# S3 method for class 'statnnet_selection'
summary(object, ...)

# S3 method for class 'statnnet_selection'
coef(object, ...)

# S3 method for class 'statnnet_selection'
predict(object, newdata, ...)

# S3 method for class 'statnnet_selection'
plot(x, ...)
```

## Arguments

- x, object:

  A `statnnet_selection` object.

- ...:

  Additional arguments passed to the underlying method where documented.

- newdata:

  Optional new predictor data. For a formula selection this is a data
  frame; for a matrix selection it may be a matrix or data frame.

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
summary(selection)
#> statnnet model-selection summary
#> Selected formula: y ~ x1
#> Architecture: 1-1-1; 4 weights
#> BIC: 81.590
#> Search: 5 candidates; 0 failed candidates; 0 failed starts
#> 
#> Terms
#>  term selected
#>    x1     TRUE
#>    x2    FALSE
coef(selection)
#>       b->h1      i1->h1        b->o       h1->o 
#>  -0.6055015  -0.2181255  15.3085551 -40.4499540 
head(predict(selection))
#>         [,1]
#> 1 -0.2601090
#> 2  1.3941326
#> 3 -0.6992796
#> 4  4.0559063
#> 5  1.6830968
#> 6 -0.6673043
plot(selection)
```
