# Selecting and interpreting a shallow neural network

The usual `statnnet` workflow selects a parsimonious network before
constructing regression-style summaries. The example uses the bundled
Auto MPG data and a small search so that the vignette can be rebuilt
quickly.

``` r

library(nnet)
library(statnnet)

data("auto_mpg", package = "statnnet")

analysis_data <- auto_mpg
continuous_variables <- c(
  "mpg", "displacement", "horsepower", "weight", "acceleration",
  "model_year"
)
analysis_data[continuous_variables] <- scale(
  analysis_data[continuous_variables]
)
```

## Select inputs and hidden-layer size

[`select_nnet()`](https://andrew-mcinerney.github.io/statnnet/reference/select_nnet.md)
alternates between hidden-layer and input selection using an information
criterion. Formula terms are selected as groups, so the indicator
columns for a factor enter or leave together. Multiple starting values
are used for every candidate model; only converged starts are eligible
for selection.

``` r

selection <- select_nnet(
  mpg ~ cylinders + displacement + horsepower + weight + acceleration +
    model_year + origin,
  data = analysis_data,
  max_size = 2,
  n_starts = 3,
  criterion = "BIC",
  maxit = 1000,
  seed = 2026,
  progress = FALSE
)

selection
#> statnnet model selection
#> Call:
#> select_nnet(x = mpg ~ cylinders + displacement + horsepower + 
#>     weight + acceleration + model_year + origin, data = analysis_data, 
#>     max_size = 2, n_starts = 3, criterion = "BIC", maxit = 1000, 
#>     seed = 2026, progress = FALSE)
#> Selected model: 5 input columns, 2 hidden units, 15 network weights
#> BIC: 332.725; 33 candidate models; 99 converged starts
#> Selected terms: horsepower, weight, model_year, origin
#> Excluded terms: cylinders, displacement, acceleration
summary(selection)
#> statnnet model-selection summary
#> Selected formula: mpg ~ horsepower + weight + model_year + origin
#> Architecture: 5-2-1; 15 weights
#> BIC: 332.725
#> Search: 33 candidates; 0 failed candidates; 0 failed starts
#> 
#> Terms
#>          term selected
#>    horsepower     TRUE
#>        weight     TRUE
#>    model_year     TRUE
#>        origin     TRUE
#>     cylinders    FALSE
#>  displacement    FALSE
#>  acceleration    FALSE
```

The evaluated search path is retained and can be inspected or plotted.

``` r

head(selection_path(selection))
#>   evaluation  phase step           candidate
#> 1          1 hidden    1              size_1
#> 2          2 hidden    1              size_2
#> 3          3  input    1    remove_cylinders
#> 4          4  input    1 remove_displacement
#> 5          5  input    1   remove_horsepower
#> 6          6  input    1       remove_weight
#>                                                                        included_terms
#> 1 cylinders + displacement + horsepower + weight + acceleration + model_year + origin
#> 2 cylinders + displacement + horsepower + weight + acceleration + model_year + origin
#> 3             displacement + horsepower + weight + acceleration + model_year + origin
#> 4                cylinders + horsepower + weight + acceleration + model_year + origin
#> 5              cylinders + displacement + weight + acceleration + model_year + origin
#> 6          cylinders + displacement + horsepower + acceleration + model_year + origin
#>   n_terms n_inputs size n_parameters criterion criterion_value converged_starts
#> 1       7       11    1           15       BIC        362.8806                3
#> 2       7       11    2           28       BIC        363.5503                3
#> 3       6        7    1           11       BIC        358.7954                3
#> 4       6       10    1           14       BIC        357.0597                3
#> 5       6       10    1           14       BIC        822.0933                3
#> 6       6       10    1           14       BIC        398.4909                3
#>   n_starts best_start selected
#> 1        3          3    FALSE
#> 2        3          3    FALSE
#> 3        3          3    FALSE
#> 4        3          3    FALSE
#> 5        3          1    FALSE
#> 6        3          3    FALSE
plot(selection)
```

![BIC values for the candidate neural-network models evaluated during
selection.](selection-workflow_files/figure-html/path-1.png)

## Check the numerical diagnostics

Converting the result with
[`statnnet()`](https://andrew-mcinerney.github.io/statnnet/reference/statnnet.md)
checks the selected `nnet` fit and reports whether covariance-based
summaries are numerically available. Selection is unpenalised, and a
converged fit may still have poorly conditioned curvature.

Numerical availability does not show that Wald p-values or delta-method
intervals have their nominal coverage. The Gaussian simulations
motivating the method favour parsimonious networks, including a hidden
layer no larger than the fitted input layer.
[`statnnet()`](https://andrew-mcinerney.github.io/statnnet/reference/statnnet.md)
records this comparison and warns when there are more hidden nodes than
input nodes. The inferential quantities remain approximate and
conditional on the setting studied.

The package can calculate the same summaries for supported
binary-response fits, but the package authors have not established their
frequentist calibration; binary Wald tests and intervals should
therefore be treated as exploratory.

``` r

selected_model <- suppressWarnings(statnnet(selection))
selected_model
#> statnnet: statistical summary of a fitted nnet model
#> Formula: mpg ~ horsepower + weight + model_year + origin
#> Architecture: 5-2-1; 15 weights; n = 392
#> Response: continuous; decay: 0
#> Selection: BIC = 332.725; uncertainty is conditional on the selected model
#> nnet convergence code: 0
#> Covariance: unavailable - the penalised information matrix is not positive definite
```

Here we refit the selected architecture with a small, pre-specified
decay to demonstrate the inferential interface. This decay value is not
selected by BIC. In an analysis it should be justified and examined in a
sensitivity analysis, while the resulting uncertainty remains
conditional on the selected structure, penalty, and local optimum.

``` r

regularised_fit <- nnet(
  selection$formula,
  data = analysis_data,
  size = selection$size,
  linout = TRUE,
  decay = 0.02,
  Hess = TRUE,
  Wts = selection$fit$wts,
  maxit = 1000,
  trace = FALSE
)

model <- statnnet(
  regularised_fit,
  data = analysis_data,
  response = "continuous"
)

model
#> statnnet: statistical summary of a fitted nnet model
#> Formula: mpg ~ horsepower + weight + model_year + origin
#> Architecture: 5-2-1; 15 weights; n = 392
#> Response: continuous; decay: 0.02
#> nnet convergence code: 0
#> Covariance: available (reciprocal condition number 2.27e-05)
summary(model, effects = "none")
#> statnnet summary
#> Formula: mpg ~ horsepower + weight + model_year + origin
#> Response: continuous; architecture: 5-2-1
#> 
#> Effect summaries: omitted
#> 
#> Grouped input-to-hidden Wald summaries
#>        term             model_columns n_weights effective_df wald_chisq
#>  horsepower                horsepower         2     1.889432   16.61118
#>      weight                    weight         2     1.969374   50.47619
#>  model_year                model_year         2     1.982435   67.15539
#>      origin originEurope, originJapan         4     3.899817   13.73133
#>       p_value
#>  2.111163e-04
#>  1.031783e-11
#>  2.521446e-15
#>  7.548034e-03
anova(model)
#> Grouped input-to-hidden Wald summaries
#>        term             model_columns n_weights effective_df wald_chisq
#>  horsepower                horsepower         2     1.889432   16.61118
#>      weight                    weight         2     1.969374   50.47619
#>  model_year                model_year         2     1.982435   67.15539
#>      origin originEurope, originJapan         4     3.899817   13.73133
#>       p_value
#>  2.111163e-04
#>  1.031783e-11
#>  2.521446e-15
#>  7.548034e-03
```

## Describe the fitted response surface

[`pce()`](https://andrew-mcinerney.github.io/statnnet/reference/pce.md)
calculates finite differences of the partial-dependence function. These
are descriptive summaries of the fitted network rather than causal
effects, and their pointwise intervals are conditional on the modelling
choices above.

``` r

pce(model, "weight", length_out = 5)
#>   variable contrast       value step     estimate  std_error   conf_low
#> 1   weight     <NA> -1.60652233    1 -0.829196198 0.08958951 -1.0047884
#> 2   weight     <NA> -0.81843976    1 -0.730285350 0.05891874 -0.8457640
#> 3   weight     <NA> -0.03035718    1 -0.478778741 0.05854461 -0.5935241
#> 4   weight     <NA>  0.75772540    1 -0.209584604 0.09348839 -0.3928185
#> 5   weight     <NA>  1.54580798    1 -0.009546099 0.11030917 -0.2257481
#>     conf_high
#> 1 -0.65360399
#> 2 -0.61480675
#> 3 -0.36403341
#> 4 -0.02635073
#> 5  0.20665590
```

``` r

plot(model, "weight", length_out = 25)
```

![Partial covariate-effect curve for standardized vehicle weight in the
selected Auto MPG
network.](selection-workflow_files/figure-html/pce-plot-1.png)

The selected model retains the original `nnet` fit, so ordinary
predictions remain available through
[`predict()`](https://rdrr.io/r/stats/predict.html).

``` r

head(predict(model, newdata = analysis_data))
#>        [,1]
#> 1 -1.035954
#> 2 -1.134669
#> 3 -1.029109
#> 4 -1.027883
#> 5 -1.024023
#> 6 -1.296471
```
