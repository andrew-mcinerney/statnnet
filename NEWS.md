# statnnet 0.1.0

* Provides `select_nnet()` for information-criterion selection of input terms and
  hidden-layer size using multiple `nnet` initialisations.
* Provides `statnnet_selection` objects with print, summary, coefficient,
  prediction, plot, and `selection_path()` methods.
* `statnnet()` accepts a model-selection result and records that subsequent
  uncertainty summaries are conditional on the selected model.
* Provides `plot_network()` for visualising fitted network weights and
  inferential annotations.
* `statnnet()` retains and validates the fitted `nnet` object.
* Reports explicit objective, Hessian, covariance, convergence, and conditioning
  diagnostics.
* Reports inferential-scope diagnostics: `statnnet()` warns when hidden nodes
  outnumber fitted input nodes and flags binary-response calibration as
  unvalidated rather than equating numerical invertibility with calibration.
* Provides structured grouped Wald summaries and partial covariate effects.
* Includes the openly licensed Auto MPG data, a reproducible data-preparation
  script, and data provenance.
