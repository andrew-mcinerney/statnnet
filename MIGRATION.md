# Migrating from `interpretnn` to `statnnet`

The package is deliberately limited to supported models fitted with
`nnet::nnet()`. It now incorporates the model-selection procedure formerly
provided by the development package `selectnn`.

The new workflow is:

```r
fit <- nnet::nnet(
  response ~ .,
  data = analysis_data,
  size = 2,
  linout = TRUE,
  decay = 0.01,
  Hess = TRUE,
  trace = FALSE
)

model <- statnnet::statnnet(
  fit,
  formula = response ~ .,
  data = analysis_data,
  response = "continuous"
)
```

Alternatively, select an unpenalised Gaussian model before constructing its
statistical summary:

```r
selection <- statnnet::select_nnet(
  response ~ .,
  data = analysis_data,
  max_size = 5,
  n_starts = 10,
  seed = 123
)

model <- statnnet::statnnet(selection)
```

The former `selectnn()` name and its low-level exported phase functions are not
retained. Use `select_nnet()` for the complete procedure and
`selection_path()` to inspect its candidate models. Adapters for Keras, Torch,
Luz, `neuralnet`, and other fitting packages are not part of `statnnet`.
