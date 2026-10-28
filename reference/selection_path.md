# Extract the model-selection path

Extract the model-selection path

## Usage

``` r
selection_path(object)
```

## Arguments

- object:

  A `statnnet_selection` object.

## Value

A data frame with one row per candidate model evaluated during the
selection procedure.

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
selection_path(selection)
#>   evaluation       phase step candidate included_terms n_terms n_inputs size
#> 1          1      hidden    1    size_1        x1 + x2       2        2    1
#> 2          2       input    1 remove_x1             x2       1        1    1
#> 3          3       input    1 remove_x2             x1       1        1    1
#> 4          4 fine_hidden    1    size_1             x1       1        1    1
#> 5          5 fine_hidden    1    size_2             x1       1        1    2
#>   n_parameters criterion criterion_value converged_starts n_starts best_start
#> 1            6       BIC        85.13239                2        2          2
#> 2            5       BIC       256.65938                2        2          2
#> 3            5       BIC        81.59030                2        2          2
#> 4            5       BIC        81.67392                2        2          1
#> 5            8       BIC        91.20925                2        2          2
#>   selected
#> 1    FALSE
#> 2    FALSE
#> 3     TRUE
#> 4    FALSE
#> 5    FALSE
```
