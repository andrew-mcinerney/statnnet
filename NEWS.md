# statnnet 0.1.0.9000

* Refactored the package as an `nnet`-only statistical companion.
* Renamed the constructor from `interpretnn()` to `statnnet()` so it matches
  the package name. This is a breaking API change.
* `statnnet()` now retains and validates an existing fitted `nnet` object.
* Added explicit objective, Hessian, covariance, convergence, and conditioning
  diagnostics.
* Added structured grouped Wald summaries and partial covariate effects.
* Removed model fitting, BIC/model selection, and other neural-network backends.
