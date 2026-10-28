#' Augment a fitted nnet model with statistical summaries
#'
#' `statnnet()` adds covariance estimates, Wald summaries, and covariate-effect
#' calculations to a supported model fitted by [nnet::nnet()]. The original
#' fit is retained. A model selected by [select_nnet()] can also be supplied.
#'
#' Covariance calculations are returned only when the fitted Hessian passes the
#' numerical checks. Wald tests and delta-method intervals are still
#' approximations conditional on the fitted architecture and local solution.
#' A warning is given when the hidden layer is larger than the fitted input
#' layer. Binary-response summaries should be treated as exploratory because
#' their frequentist calibration has not been studied by the package authors.
#'
#' @param object A fitted object inheriting from `"nnet"`, or a
#'   `"statnnet_selection"` object returned by [select_nnet()].
#' @param formula The formula used to fit `object`. This may be omitted for an
#'   `nnet.formula` object.
#' @param data The original training data.
#' @param response Either `"continuous"` or `"binary"`. When `NULL`, the
#'   response type is inferred from the fitted `nnet` object.
#' @param covariance_tol Relative numerical tolerance used to decide whether
#'   the penalised information matrix is reliably invertible.
#' @param objective_tol Relative tolerance for validating the reconstructed
#'   fitting criterion against `object$value`.
#' @param ... Reserved for methods.
#'
#' @return An object of class `"statnnet"` containing the original fit,
#'   reconstructed model information, covariance calculations, and diagnostics.
#' @references
#' McInerney, A. and Burke, K. (2026). Investigating Statistical Inference and
#' Covariate Effects in Shallow Neural Networks. arXiv:2311.08139.
#' \doi{10.48550/arXiv.2311.08139}.
#' @examples
#' data("auto_mpg", package = "statnnet")
#' example_data <- auto_mpg
#' example_data$mpg <- as.numeric(scale(example_data$mpg))
#' example_data$weight <- as.numeric(scale(example_data$weight))
#' set.seed(3)
#' fit <- nnet::nnet(
#'   mpg ~ weight + origin,
#'   data = example_data,
#'   size = 1,
#'   linout = TRUE,
#'   decay = 0.1,
#'   Hess = TRUE,
#'   maxit = 500,
#'   trace = FALSE
#' )
#' model <- statnnet(fit, data = example_data, response = "continuous")
#' model
#' @export
statnnet <- function(object, ...) {
  UseMethod("statnnet")
}

#' @rdname statnnet
#' @export
statnnet.default <- function(object, ...) {
  stop(
    paste0(
      "`object` must inherit from \"nnet\" or \"statnnet_selection\"; ",
      "statnnet supports no other fitting backends."
    ),
    call. = FALSE
  )
}

#' @rdname statnnet
#' @export
statnnet.statnnet_selection <- function(object, covariance_tol = 1e-10,
                                        objective_tol = 1e-6, ...) {
  if (length(list(...)) > 0L) {
    stop("Unused arguments were supplied to `statnnet()`.", call. = FALSE)
  }
  result <- statnnet.nnet(
    object = object$fit,
    formula = object$formula,
    data = object$training_data,
    response = "continuous",
    covariance_tol = covariance_tol,
    objective_tol = objective_tol
  )
  call <- match.call()
  call[[1L]] <- quote(statnnet)
  result$call <- call
  result$selection <- list(
    criterion = object$criterion,
    criterion_value = object$criterion_value,
    selected_terms = object$selected_terms,
    excluded_terms = object$excluded_terms,
    size = object$size,
    search_path = object$search_path,
    diagnostics = object$diagnostics,
    call = object$call
  )
  result
}

#' @rdname statnnet
#' @export
statnnet.nnet <- function(object, formula = NULL, data = NULL,
                          response = NULL, covariance_tol = 1e-10,
                          objective_tol = 1e-6, ...) {
  if (!inherits(object, "nnet")) {
    stop("`object` must inherit from \"nnet\".", call. = FALSE)
  }
  if (length(list(...)) > 0L) {
    stop("Unused arguments were supplied to `statnnet()`.", call. = FALSE)
  }
  if (is.null(data)) {
    stop("`data` must contain the original training data.", call. = FALSE)
  }
  if (!is.numeric(covariance_tol) || length(covariance_tol) != 1L ||
      !is.finite(covariance_tol) || covariance_tol <= 0) {
    stop("`covariance_tol` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(objective_tol) || length(objective_tol) != 1L ||
      !is.finite(objective_tol) || objective_tol <= 0) {
    stop("`objective_tol` must be one positive finite number.", call. = FALSE)
  }

  .validate_nnet_structure(object)

  if (is.null(formula)) {
    if (is.null(object$terms)) {
      stop(
        "`formula` is required because it cannot be recovered from this nnet object.",
        call. = FALSE
      )
    }
    formula <- stats::formula(object$terms)
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }

  inferred_response <- .nnet_response_type(object)
  if (is.null(response)) {
    response <- inferred_response
  } else {
    response <- match.arg(response, c("continuous", "binary"))
    if (!identical(response, inferred_response)) {
      stop(
        sprintf(
          "`response = \"%s\"` is incompatible with the fitted nnet output criterion.",
          response
        ),
        call. = FALSE
      )
    }
  }

  # Reconstruct the data used by nnet before calculating any summaries.
  reconstructed <- .reconstruct_nnet_data(
    object = object,
    formula = formula,
    data = data,
    response = response,
    evaluation_environment = parent.frame()
  )

  # Check that the reconstructed design matrix gives the fitted predictions.
  weights <- as.numeric(object$wts)
  fitted_internal <- .nn_predict_matrix(
    reconstructed$x,
    weights,
    hidden = object$n[2L],
    response = response
  )
  fitted_nnet <- as.numeric(object$fitted.values)
  prediction_error <- max(abs(fitted_internal - fitted_nnet))
  prediction_scale <- max(1, max(abs(fitted_nnet)))
  if (!is.finite(prediction_error) ||
      prediction_error > objective_tol * prediction_scale) {
    stop(
      "The reconstructed model matrix does not reproduce the fitted nnet predictions.",
      call. = FALSE
    )
  }

  # Reconstruct the fitting criterion and its curvature.
  objective <- .validate_nnet_objective(
    object = object,
    y = reconstructed$y,
    fitted = fitted_internal,
    response = response,
    tolerance = objective_tol
  )
  curvature <- .nnet_curvature(
    object = object,
    x = reconstructed$x,
    y = reconstructed$y,
    response = response,
    rss = objective$unpenalised,
    tolerance = covariance_tol
  )

  # Use nnet's coefficient labels when they are available.
  parameter_names <- names(stats::coef(object))
  if (is.null(parameter_names) || length(parameter_names) != length(weights)) {
    parameter_names <- paste0("weight_", seq_along(weights))
  }
  names(weights) <- parameter_names
  dimnames(curvature$hessian_penalised) <- list(parameter_names, parameter_names)
  dimnames(curvature$hessian_unpenalised) <- list(parameter_names, parameter_names)
  dimnames(curvature$information_penalised) <- list(parameter_names, parameter_names)
  dimnames(curvature$information_unpenalised) <- list(parameter_names, parameter_names)
  if (!is.null(curvature$covariance)) {
    dimnames(curvature$covariance) <- list(parameter_names, parameter_names)
    dimnames(curvature$shrinkage) <- list(parameter_names, parameter_names)
  }

  groups <- .make_term_groups(
    reconstructed$assign,
    reconstructed$term_labels,
    reconstructed$column_names,
    n_inputs = object$n[1L],
    n_hidden = object$n[2L]
  )

  # Collect numerical and model-architecture diagnostics in one place.
  diagnostics <- curvature$diagnostics
  diagnostics$convergence <- object$convergence
  diagnostics$objective_agrees <- objective$agrees
  diagnostics$objective_difference <- objective$difference
  diagnostics$prediction_max_abs_difference <- prediction_error
  diagnostics$architecture <- unname(as.integer(object$n))
  diagnostics$n_parameters <- length(weights)
  diagnostics$n_observations <- nrow(reconstructed$x)
  diagnostics$decay <- as.numeric(object$decay)
  diagnostics$n_input_nodes <- as.integer(object$n[1L])
  diagnostics$n_hidden_nodes <- as.integer(object$n[2L])
  diagnostics$hidden_to_input_ratio <- (
    diagnostics$n_hidden_nodes / diagnostics$n_input_nodes
  )
  diagnostics$architecture_calibration_warning <- (
    diagnostics$n_hidden_nodes > diagnostics$n_input_nodes
  )
  diagnostics$binary_calibration_unvalidated <- identical(response, "binary")

  if (!identical(as.integer(object$convergence), 0L)) {
    diagnostics$covariance_available <- FALSE
    diagnostics$covariance_reason <- paste(
      "the nnet optimiser did not converge (code",
      as.integer(object$convergence),
      ")"
    )
    curvature$covariance <- NULL
    curvature$shrinkage <- NULL
  }

  # Keep the original fit so the usual nnet methods remain available.
  result <- list(
    fit = object,
    call = match.call(),
    formula = formula,
    terms = reconstructed$terms,
    contrasts = reconstructed$contrasts,
    xlevels = reconstructed$xlevels,
    training_data = reconstructed$data,
    x = reconstructed$x,
    y = reconstructed$y,
    response = response,
    response_levels = reconstructed$response_levels,
    weights = weights,
    decay = as.numeric(object$decay),
    architecture = list(
      inputs = as.integer(object$n[1L]),
      hidden = as.integer(object$n[2L]),
      outputs = as.integer(object$n[3L])
    ),
    objective = objective,
    hessian_penalised = curvature$hessian_penalised,
    hessian_unpenalised = curvature$hessian_unpenalised,
    information_penalised = curvature$information_penalised,
    information_unpenalised = curvature$information_unpenalised,
    covariance = curvature$covariance,
    shrinkage = curvature$shrinkage,
    groups = groups,
    variables = reconstructed$variables,
    diagnostics = diagnostics,
    covariance_tol = covariance_tol,
    selection = NULL
  )
  class(result) <- "statnnet"

  # Warn once, after all diagnostic checks have been completed.
  diagnostic_warnings <- character()
  if (!isTRUE(diagnostics$covariance_available)) {
    diagnostic_warnings <- c(
      diagnostic_warnings,
      paste0(
        "Covariance-based summaries are unavailable: ",
        diagnostics$covariance_reason,
        "."
      )
    )
  } else if (isTRUE(diagnostics$conditioning_warning)) {
    diagnostic_warnings <- c(
      diagnostic_warnings,
      sprintf(
        "The penalised information matrix is poorly conditioned (reciprocal condition number %.3g).",
        diagnostics$information_rcond
      )
    )
  }
  if (isTRUE(diagnostics$architecture_calibration_warning)) {
    diagnostic_warnings <- c(
      diagnostic_warnings,
      sprintf(
        paste(
          "The network has more hidden nodes (%d) than fitted input nodes (%d).",
          "Numerical covariance availability does not establish calibrated",
          "inference; use a more parsimonious architecture or treat Wald and",
          "delta-method results as exploratory."
        ),
        diagnostics$n_hidden_nodes,
        diagnostics$n_input_nodes
      )
    )
  }
  if (length(diagnostic_warnings)) {
    warning(paste(diagnostic_warnings, collapse = " "), call. = FALSE)
  }

  result
}
