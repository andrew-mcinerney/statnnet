# statnnet 0.1.0.9000

* Refactored the package as an `nnet`-only statistical companion.
* `interpretnn()` now retains and validates an existing fitted `nnet` object.
* Added explicit objective, Hessian, covariance, convergence, and conditioning
  diagnostics.
* Added structured grouped Wald summaries and partial covariate effects.
* Removed model fitting, BIC/model selection, and other neural-network backends.
