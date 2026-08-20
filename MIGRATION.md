# Migrating from `interpretnn` to `statnnet`

The package is now deliberately limited to supported models fitted with
`nnet::nnet()`. Both the package and its constructor are called `statnnet`,
and the returned object has class `"statnnet"`.

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

The former fitting and model-selection interface (`nn_fit()` and BIC fields),
and adapters for Keras, Torch, Luz, `neuralnet`, and other packages, are not
part of `statnnet`. Fit and validate the desired `nnet` model first, then pass
that unchanged fit to `statnnet()`.
