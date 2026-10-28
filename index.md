# statnnet

`statnnet` provides model selection and regression-style summaries for
small, single-hidden-layer neural networks fitted with `nnet`. It can
select an unpenalised Gaussian model by information criterion or work
from an existing supported `nnet` fit.

## Installation

``` r

# install.packages("pak")
pak::pak("andrew-mcinerney/statnnet")
```

## Basic workflow

The main workflow chooses the input terms and hidden-layer size from
multiple random starts, then adds diagnostics and regression-style
summaries:

``` r

library(nnet)
library(statnnet)

data("auto_mpg", package = "statnnet")
analysis_data <- auto_mpg
continuous_variables <- c(
  "mpg", "displacement", "horsepower", "weight",
  "acceleration", "model_year"
)
analysis_data[continuous_variables] <- lapply(
  analysis_data[continuous_variables],
  function(x) as.numeric(scale(x))
)
```

``` r

selection <- select_nnet(
  mpg ~ displacement + horsepower + weight + acceleration +
    model_year + origin,
  data = analysis_data,
  max_size = 3,
  n_starts = 10,
  seed = 2026
)

selection
selection_path(selection)
model <- statnnet(selection)
```

An existing supported fit can also be augmented directly:

``` r

set.seed(2026)
fit <- nnet(
  mpg ~ weight + model_year + origin,
  data = analysis_data,
  size = 1,
  linout = TRUE,
  decay = 0.1,
  Hess = TRUE,
  maxit = 1000,
  trace = FALSE
)

model <- statnnet(
  fit,
  data = analysis_data,
  response = "continuous"
)

model
#> statnnet: statistical summary of a fitted nnet model
#> Formula: mpg ~ weight + model_year + origin
#> Architecture: 4-1-1; 7 weights; n = 392
#> Response: continuous; decay: 0.1
#> nnet convergence code: 0
#> Covariance: available (reciprocal condition number 0.000405)
anova(model)
#> Grouped input-to-hidden Wald summaries
#>        term             model_columns n_weights effective_df wald_chisq
#>      weight                    weight         1    0.9926738  121.79808
#>  model_year                model_year         1    0.9982324  110.91756
#>      origin originEurope, originJapan         2    1.9936731   15.20982
#>       p_value
#>  2.499232e-28
#>  6.135708e-26
#>  4.937229e-04
```

`summary(model)` returns structured average observed-row predictive
effects and grouped Wald-type summaries.
[`pce()`](https://andrew-mcinerney.github.io/statnnet/reference/pce.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) provide
finite-difference partial covariate effects with delta-method or
simulation uncertainty.

``` r

summary(model)
pce(model, "weight")
plot(model, "weight")
plot_network(model, annotate = "p_value")
```

These summaries describe the fitted response function. They are
conditional on the chosen network and are not post-selection or causal
inference.

## Citation

Run `citation("statnnet")` for formatted references and BibTeX entries.
If you use
[`select_nnet()`](https://andrew-mcinerney.github.io/statnnet/reference/select_nnet.md)
for input or hidden-layer selection, cite McInerney and Burke (2025). If
you use covariance estimates, Wald summaries, grouped tests, or
covariate-effect methods, cite McInerney and Burke (2026). Analyses
using both parts of the package should cite both papers.

## Inference and diagnostics

When a covariance matrix is available, the fitted Hessian has passed the
package’s numerical checks. The method is intended for small,
parsimonious networks.
[`statnnet()`](https://andrew-mcinerney.github.io/statnnet/reference/statnnet.md)
records the hidden-to-input size ratio and reports an architecture
diagnostic when the hidden layer is larger than the fitted input layer.

Wald p-values and delta-method intervals are approximate and conditional
on the architecture, penalty, and local optimum.
