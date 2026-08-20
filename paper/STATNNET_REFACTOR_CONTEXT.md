# `statnnet` Refactor Context

## Purpose of this document

This document is the handoff specification for refactoring the existing `interpretnn` R package into `statnnet`.
It records the intended package scope, statistical conventions, public workflow, paper claims, migration sequence, and completion criteria.

The corresponding R Journal draft is `new-draft.Rmd`, titled:

> `statnnet`: Regression-Style Statistical Summaries for Neural Networks Fitted with `nnet`

Treat the paper as the specification for the intended user-facing behaviour, but do not assume that the current package already implements every claim in the draft.
Inspect the package before editing, identify differences between the code and this document, and report any statistical or API ambiguity rather than silently choosing a different method.

Do not edit generated `.tex`, `.pdf`, `.html`, or documentation files when their source should be edited instead.

## Authoritative references

Use the following sources, in order, when checking alignment:

1. The intended scope and user-facing workflow in `new-draft.Rmd`.
2. The installed `nnet` documentation and source for the fitting criterion, parameter ordering, weight decay, Hessian, and prediction behaviour.
3. McInerney and Burke, *Feedforward neural networks as statistical models: Improving interpretability through uncertainty quantification*, for the covariance, Wald, and PCE methodology.
4. The refactored package tests, which should encode verified behaviour rather than preserve accidental legacy behaviour.

When these sources disagree, stop and document the discrepancy before changing statistically material code.

## Refactor objective

Refactor `interpretnn` into a deliberately narrow statistical companion to the `nnet` package.

The package should:

- accept supported fitted `nnet` objects;
- retain the original fit rather than introduce another optimiser;
- reconstruct and validate the statistical model associated with the fit;
- provide covariance estimates and regression-style summaries;
- provide individual-parameter and grouped covariate-level Wald-type summaries;
- provide partial covariate-effect (PCE) estimates and plots with uncertainty;
- expose convergence and covariance diagnostics;
- use familiar S3 methods where appropriate;
- preserve prediction by delegating to the underlying `nnet` fit.

The central package identity is:

> Statistical summaries for small, parsimonious, single-hidden-layer neural networks fitted using `nnet`.

## Explicit non-goals

The refactor must not present `statnnet` as:

- a general neural-network inference framework;
- an interface to `neuralnet`, `keras`, `torch`, `luz`, or other fitting backends;
- a package for deep, convolutional, recurrent, transformer, or highly overparameterised models;
- a new neural-network optimiser;
- a model-selection package;
- a BIC-based package or a continuation of earlier model-selection functionality;
- a provider of post-selection inference;
- a causal-inference package;
- a replacement for predictive validation or model checking.

Do not add backend abstractions or speculative extensibility for deep-learning frameworks during this refactor.

## Intended public workflow

The current paper retains `interpretnn()` as the constructor inside the renamed `statnnet` package.
Do not mechanically rename this function to `statnnet()` unless the paper and package API are deliberately changed together.

The intended workflow is:

```r
library(nnet)
library(statnnet)

nn <- nnet(
  response ~ .,
  data = analysis_data,
  size = 2,
  linout = TRUE,
  decay = 0.01,
  Hess = TRUE,
  trace = FALSE
)

snn <- interpretnn(
  object = nn,
  formula = response ~ .,
  data = analysis_data,
  response = "continuous"
)

print(snn)
summary(snn)
anova(snn)
vcov(snn)
plot(snn)
plotnn(snn)
predict(snn, newdata = new_data)
```

The intended augmented class is `"statnnet"`.
The original fitted `nnet` object should be retained internally.

If backward compatibility is required, provide a documented migration route from the old package and class names.
Do not preserve old names through hidden or ambiguous dispatch behaviour.

## Supported-model contract

Initially support only models satisfying all of the following conditions:

- The object inherits from `"nnet"` or an explicitly supported `nnet` formula subclass.
- There is exactly one hidden layer.
- There is exactly one output node.
- There are no skip-layer connections.
- The fit is either:
  - continuous-response regression with a linear output; or
  - binary-response regression with a single logistic output and entropy fitting.
- The response, model matrix, fitted architecture, and stored weights have compatible dimensions.
- The original formula and training data are available when they cannot be reconstructed unambiguously from the fitted object.
- The scalar `decay` convention described below is used.

Initially reject, with informative messages:

- softmax models;
- censored outputs;
- multinomial responses;
- multiple output nodes;
- skip-layer connections;
- unsupported masks or fixed-weight configurations;
- objects from other neural-network packages;
- non-finite weights, responses, design matrices, or Hessians.

The present paper mathematics is unweighted.
Either reject non-uniform case weights initially or extend the implementation, tests, and paper equations together before claiming support.

## Formula, factors, and model-matrix handling

Factor handling is part of the statistical contribution and must be explicit.

- Reconstruct the design matrix using the original formula, data, contrasts, and factor levels.
- Preserve the mapping from each design-matrix column to its original covariate.
- Use that mapping for grouped factor tests in `anova()`.
- Store or reconstruct the information needed to apply identical transformations to `newdata`.
- Check that the reconstructed design matrix agrees with the number and order of input weights in the fitted `nnet` object.
- Do not infer factor groupings only from similar column names.
- Add tests for reordered factor levels, absent levels in `newdata`, and non-default contrasts.

The paper currently describes preprocessed insurance data while also discussing grouped factor tests.
The final implementation and example must resolve this by either retaining factors until model-matrix construction or explicitly supplying grouping metadata.

## Exact `nnet` fitting criteria

Use notation consistent with the paper.
Let:

- `theta` be the vector of all `K` fitted `nnet` weights;
- `delta >= 0` be the scalar `decay` value;
- `NN(x_i, theta)` be the network output.

For scalar `decay`, `nnet` penalises all fitted weights, including hidden- and output-layer intercept weights.

### Continuous response

For a linear output, `nnet` minimises:

```text
Q_delta,G(theta)
  = sum_i {y_i - NN(x_i, theta)}^2
    + delta * theta' theta.
```

The stored `value` must equal residual sum of squares plus the weight-decay term, up to numerical tolerance.

Estimate the Gaussian residual variance using the convention in the paper:

```text
sigma_hat^2 = RSS / n.
```

This is the maximum-likelihood divisor.
If a residual-degrees-of-freedom divisor is desired instead, change the paper, implementation, and tests together.

### Binary response

For a binary factor response, `nnet` uses one logistic output with entropy fitting and minimises:

```text
Q_delta,B(theta)
  = -sum_i [y_i log{NN(x_i, theta)}
            + (1 - y_i) log{1 - NN(x_i, theta)}]
    + delta * theta' theta.
```

The stored `value` must equal negative Bernoulli log-likelihood plus the weight-decay term, up to numerical tolerance.

## Hessian and covariance calculation

When `Hess = TRUE`, the `nnet` Hessian is the Hessian of the penalised fitting criterion.
Write:

```text
H_delta = d^2 Q_delta(theta_hat) / d theta d theta'.
```

For scalar `decay`, recover the unpenalised curvature on the native `nnet` scale using:

```text
H_0 = H_delta - 2 * delta * I_K.
```

Convert to the likelihood scale as follows.

For binary entropy fitting:

```text
J_0     = H_0
J_delta = H_delta.
```

For Gaussian regression:

```text
J_0     = H_0     / (2 * sigma_hat^2)
J_delta = H_delta / (2 * sigma_hat^2).
```

Estimate the covariance matrix of the penalised estimate using:

```text
Cov_hat(theta_hat)
  = inverse(J_delta) %*% J_0 %*% inverse(J_delta).
```

When `decay = 0`, this reduces to the inverse observed-information approximation.

Do not silently replace an unstable inverse with an arbitrary generalised inverse, ridge adjustment, or nearest-positive-definite matrix.
If an alternative is added, expose it as an explicit method, document the statistical interpretation, and test it separately.

## Wald-type summaries

### Individual weights

For parameter `theta_m`, compute:

```text
W_m = theta_hat_m^2 / Cov_hat(theta_hat)[m, m]
```

and compare it with a chi-squared reference distribution with one degree of freedom.
Present individual-weight results primarily as diagnostics rather than as the main scientific interpretation.

### Input-level and covariate-level groups

For input `j`, let `omega_j` contain its `q` input-to-hidden weights and let `S_j` select these weights from `theta`.
Compute:

```text
Sigma_omega_j = S_j %*% Cov_hat(theta_hat) %*% t(S_j)

W_j = t(omega_hat_j) %*%
      inverse(Sigma_omega_j) %*%
      omega_hat_j.
```

Define the shrinkage matrix:

```text
A = inverse(J_delta) %*% J_0.
```

Use effective degrees of freedom:

```text
q_tilde_j = trace(S_j %*% A %*% t(S_j)).
```

For a factor represented by `l` dummy variables, use one selection matrix for all `l * q` associated input-to-hidden weights.
The effective degrees of freedom are computed using the same trace expression.

## Partial covariate effects

Implement the PCE as the finite difference of the partial dependence function.
For a `d`-unit increase in covariate `j`:

```text
PCE_j(x, d)
  = mean_i NN(x_i with x_j = x + d, theta_hat)
    - mean_i NN(x_i with x_j = x, theta_hat).
```

Support:

- continuous covariates over an observed range;
- binary covariates as the change from 0 to 1;
- a point summary obtained by averaging the PCE over observed covariate values;
- delta-method standard errors;
- simulation from the approximate sampling distribution of `theta_hat`;
- optional exploratory plots for interactions.

PCEs are summaries of the fitted response function, not causal effects.
Warn or document that correlated predictors can lead to evaluation at poorly supported covariate combinations.

## Diagnostics contract

Diagnostics are not decorative output; they determine whether covariance-based summaries are reportable.

At minimum, check and retain:

- the `nnet` convergence code;
- whether the Hessian exists or can be recomputed from the original data;
- finiteness and symmetry of the Hessian;
- eigenvalues or another documented measure of definiteness;
- reciprocal condition number or another documented conditioning diagnostic;
- whether `J_delta` can be inverted reliably;
- whether group covariance submatrices can be inverted reliably;
- the fitted architecture, number of parameters, sample size, and `decay` value.

Use clear conditions:

- error for unsupported model structures or irrecoverable model information;
- error or unavailable result for failed covariance calculations;
- warning for numerically questionable but computable cases;
- message only for non-problematic user information.

Do not report p-values or uncertainty intervals when their required covariance calculation has failed.

## S3 interface

The paper expects the following public functions and methods:

- `interpretnn()` validates and augments the fitted `nnet` object;
- `print()` gives a compact description of the fit and diagnostics;
- `summary()` gives PCE point estimates, standard errors, and grouped Wald-type summaries;
- `anova()` gives covariate-level grouped Wald-type summaries, including factor groups;
- `vcov()` returns the covariance estimate only when diagnostics are satisfied;
- `plot()` produces PCE plots and optional uncertainty intervals;
- `plotnn()` visualises the network with statistical annotations;
- `predict()` delegates to the retained `nnet` object.

Return structured objects from `summary()` and `anova()`, with separate print methods where appropriate.
Do not make users parse printed text to recover results.

## Package and repository rename

Perform the rename deliberately rather than as a blind text replacement.

Inspect and update, as applicable:

- `DESCRIPTION`, including `Package`, title, description, URLs, and bug-report URL;
- package directory and RStudio project name;
- namespace and S3 registrations;
- roxygen package documentation;
- package-level help aliases;
- test fixtures and snapshot paths;
- vignettes and README files;
- pkgdown configuration and site links;
- GitHub Actions and other CI files;
- badges and installation instructions;
- `data()` documentation and package-qualified examples;
- citation file and metadata;
- reverse references to `interpretnn` as a package name.

Do not replace the `interpretnn()` constructor name merely because the old package was called `interpretnn`.
The current paper intentionally distinguishes the package name `statnnet` from the constructor `interpretnn()`.

## Removing BIC and model-selection scope

The current paper has removed BIC entirely.
The refactored public API, documentation, examples, and vignettes should not present BIC or neural-network model selection as part of `statnnet`.

During the refactor:

- identify all BIC and selection-related functions, fields, tests, and documentation;
- determine whether they are unused internals, legacy exports, or part of another package;
- remove them from the `statnnet` public workflow;
- avoid leaving undocumented BIC fields in the main object;
- deprecate rather than abruptly remove user-facing legacy functions when compatibility is required;
- keep any separate model-selection work outside the scope of the R Journal paper.

## Refactor sequence

### Phase 1: Inventory and baseline

1. Read `DESCRIPTION`, `NAMESPACE`, all exported R functions, S3 registrations, tests, vignettes, and package-level documentation.
2. Run the existing test suite and `R CMD check` before editing.
3. Record current failures, warnings, and undocumented behaviour.
4. Identify code for non-`nnet` backends, model selection, BIC, fitting wrappers, covariance estimation, Wald tests, PCEs, and plotting.
5. Add baseline tests around behaviour that must be preserved.

### Phase 2: Establish the `statnnet` package identity

1. Rename package metadata and documentation.
2. Update imports and remove backend-specific dependencies that are no longer needed.
3. Keep dependencies minimal; explain any new dependency before adding it.
4. Preserve the `interpretnn()` constructor unless a coordinated API decision changes it.
5. Add a short migration document for existing users.

### Phase 3: Isolate and validate the `nnet` adapter

1. Create one clear internal path for extracting weights, architecture, connection ordering, output type, convergence, `decay`, and Hessian information from `nnet` objects.
2. Remove or quarantine adapters for other fitting frameworks.
3. Validate the supported-model contract before computing statistical summaries.
4. Reconstruct the formula, response, model matrix, contrasts, and factor groupings reproducibly.

### Phase 4: Refactor the statistical core

1. Implement the exact continuous and binary `nnet` objectives.
2. Implement native-scale and likelihood-scale Hessian conversions.
3. Implement the penalised sandwich covariance calculation.
4. Implement numerical diagnostics without silent statistical substitutions.
5. Implement individual and grouped Wald-type summaries with effective degrees of freedom.
6. Implement PCE point summaries, curves, and uncertainty methods.

### Phase 5: Rebuild the S3 interface

1. Define a stable `"statnnet"` object with the original `nnet` fit retained.
2. Keep internal storage separate from the printed user interface.
3. Implement and test `print`, `summary`, `anova`, `vcov`, `plot`, and `predict` methods.
4. Ensure `predict.statnnet()` agrees exactly with `predict.nnet()` for supported inputs.
5. Ensure methods return structured, testable results.

### Phase 6: Documentation and paper examples

1. Rewrite the README and main vignette around the direct `nnet` to `statnnet` workflow.
2. Show native `summary(nn)` followed by `summary(snn)`.
3. Demonstrate a grouped factor test.
4. Demonstrate a PCE with uncertainty.
5. Demonstrate a diagnostic failure or warning.
6. Ensure the code and terminology match `new-draft.Rmd` exactly.
7. Do not manually type outputs into the paper; generate them from evaluated chunks.

### Phase 7: Release checks

1. Run documentation, tests, examples, vignettes, and `R CMD check` in a clean session.
2. Run checks on the minimum supported R version where feasible.
3. Check package installation from a source tarball.
4. Render the R Journal article in both HTML and PDF formats.
5. Run the R Journal article checks.
6. Review all warnings rather than suppressing them globally.

## Required tests

Use `testthat` and fixed seeds for stochastic fits.
Keep tests small and deterministic.

### Objective tests

- For Gaussian regression, verify `object$value` equals RSS plus `decay * sum(weights^2)`.
- For binary entropy fitting, verify `object$value` equals negative Bernoulli log-likelihood plus `decay * sum(weights^2)`.
- Test both `decay = 0` and positive scalar `decay`.
- Verify that intercept weights are included in the scalar-decay term.

### Hessian and covariance tests

- Compare the `nnet` Hessian with a numerical Hessian of the native objective on a very small model.
- Verify `H_0 = H_delta - 2 * delta * I_K`.
- Verify the Gaussian likelihood-scale conversion.
- Verify the binary likelihood-scale conversion.
- Verify that the covariance is symmetric to numerical tolerance.
- Verify the zero-decay reduction to inverse observed information.
- Test near-singular and singular examples and their conditions.

### Wald tests

- Test parameter-selection matrices against the documented `nnet` weight ordering.
- Test one continuous input, one binary input, and one multi-level factor.
- Verify grouped statistics and effective degrees of freedom against direct matrix calculations.
- Verify behaviour when a group covariance matrix is singular.

### PCE tests

- Recover a constant effect in a simple linear-like fitted function.
- Test binary 0-to-1 effects.
- Compare analytic gradients with numerical derivatives.
- Compare delta-method standard errors with direct matrix calculations.
- Test simulation-based intervals with a fixed seed.
- Verify factor and `newdata` handling.

### Interface tests

- `print()` is concise and stable.
- `summary()` and `anova()` return structured objects.
- `vcov()` returns the documented covariance or an informative failure.
- `predict.statnnet()` agrees with `predict.nnet()`.
- unsupported architectures fail early with specific messages.
- non-converged fits and unavailable Hessians are handled explicitly.

## Paper-alignment checklist

Before changing a public function or statistical calculation, check the corresponding statement in `new-draft.Rmd`.

The package and paper must agree on:

- package name: `statnnet`;
- constructor: `interpretnn()`;
- augmented class: `"statnnet"`;
- supported responses and architectures;
- scalar `decay` convention;
- treatment of intercept weights;
- Gaussian residual-variance convention;
- Hessian scaling and penalty removal;
- covariance formula;
- effective degrees of freedom;
- factor grouping;
- diagnostic behaviour;
- method names and returned quantities;
- prediction delegation;
- absence of BIC and model-selection claims;
- conditional, non-post-selection interpretation of the Wald summaries;
- non-causal interpretation of PCEs.

If the implementation reveals that a paper claim is wrong or infeasible, do not distort the code to preserve the sentence.
Document the discrepancy, determine the statistically correct behaviour, and update the package and paper together.

## Commands for verification

From the package project:

```r
devtools::document()
devtools::test()
devtools::check()
```

Also build and check the source package:

```r
devtools::build()
rcmdcheck::rcmdcheck(args = "--as-cran")
```

From the paper project, once `statnnet` is installed:

```r
rmarkdown::render("new-draft.Rmd", output_format = "all")
rjtools::initial_check_article(".")
```

## Definition of done

The refactor is complete only when all of the following are true:

- The installed package is named `statnnet`.
- The main workflow begins with a supported fitted `nnet` object.
- `interpretnn()` creates a valid `"statnnet"` object without refitting the network.
- Unsupported architectures and backends are rejected clearly.
- The continuous and binary objectives agree numerically with `nnet`.
- The Hessian, covariance, Wald, effective-degrees-of-freedom, and PCE calculations are covered by tests.
- Diagnostics prevent invalid covariance-based output from being reported.
- Prediction agrees with the retained `nnet` fit.
- Factor groupings are reconstructed from explicit model information.
- BIC and model-selection functionality are absent from the main API and paper workflow.
- Package documentation and examples use the `statnnet` name consistently.
- The full test suite passes.
- `R CMD check --as-cran` has no unexplained errors, warnings, or notes.
- The paper runs from evaluated code and renders in HTML and PDF.
- The package output shown in the paper is generated rather than manually entered.
- The final paper and package use the same terminology, formulas, supported-model contract, and public methods.

## Working style for Codex

- Inspect relevant files before editing.
- Make incremental changes rather than a large blind rewrite.
- Preserve user changes and unrelated work.
- Prefer correctness and reproducibility over clever abstractions.
- Use snake case for new internal functions and variables.
- Keep functions small and testable.
- Avoid unnecessary dependencies.
- Use `set.seed()` in tests and stochastic examples.
- Check convergence, warnings, gradients, Hessians, and numerical stability.
- Do not alter the statistical method merely to make tests pass.
- Explain every statistically material change.
- After each phase, summarise what changed and provide exact verification commands.
