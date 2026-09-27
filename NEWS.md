# statnnet 0.1.0

* Added `select_nnet()` for information-criterion selection of input terms and
  hidden-layer size using multiple `nnet` initialisations.
* Added `statnnet_selection` objects with print, summary, coefficient,
  prediction, plot, and `selection_path()` methods.
* `statnnet()` now accepts a model-selection result and records that subsequent
  uncertainty summaries are conditional on the selected model.
* Renamed `plotnn()` to `plot_network()`.
* Refactored the package as an `nnet`-only statistical companion.
* Renamed the constructor from `interpretnn()` to `statnnet()` so it matches
  the package name. This is a breaking API change.
* `statnnet()` now retains and validates an existing fitted `nnet` object.
* Added explicit objective, Hessian, covariance, convergence, and conditioning
  diagnostics.
* Added inferential-scope diagnostics: `statnnet()` now warns when hidden nodes
  outnumber fitted input nodes and flags binary-response calibration as
  unvalidated rather than equating numerical invertibility with calibration.
* Added structured grouped Wald summaries and partial covariate effects.
* Replaced the legacy Boston data and examples with the openly licensed Auto
  MPG data, including a reproducible data-preparation script and provenance.
* Removed other neural-network backends.
