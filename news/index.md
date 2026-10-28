# Changelog

## statnnet 0.1.0

- Added
  [`select_nnet()`](https://andrew-mcinerney.github.io/statnnet/reference/select_nnet.md)
  for information-criterion selection of input terms and hidden-layer
  size using multiple `nnet` initialisations.
- Added `statnnet_selection` objects with print, summary, coefficient,
  prediction, plot, and
  [`selection_path()`](https://andrew-mcinerney.github.io/statnnet/reference/selection_path.md)
  methods.
- [`statnnet()`](https://andrew-mcinerney.github.io/statnnet/reference/statnnet.md)
  now accepts a model-selection result and records that subsequent
  uncertainty summaries are conditional on the selected model.
- Added
  [`plot_network()`](https://andrew-mcinerney.github.io/statnnet/reference/plot_network.md)
  for visualising fitted network weights and inferential annotations.
- [`statnnet()`](https://andrew-mcinerney.github.io/statnnet/reference/statnnet.md)
  now retains and validates an existing fitted `nnet` object.
- Added explicit objective, Hessian, covariance, convergence, and
  conditioning diagnostics.
- Added inferential-scope diagnostics:
  [`statnnet()`](https://andrew-mcinerney.github.io/statnnet/reference/statnnet.md)
  now warns when hidden nodes outnumber fitted input nodes and flags
  binary-response calibration as unvalidated rather than equating
  numerical invertibility with calibration.
- Added structured grouped Wald summaries and partial covariate effects.
- Replaced the legacy Boston data and examples with the openly licensed
  Auto MPG data, including a reproducible data-preparation script and
  provenance.
