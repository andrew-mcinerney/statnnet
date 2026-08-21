#' Methods for statnnet objects
#'
#' @param x,object A `statnnet` object.
#' @param newdata Optional data passed unchanged to the retained `nnet` fit.
#' @param ... Additional arguments passed to the underlying method where
#'   documented.
#' @param weights Include individual-weight Wald diagnostics in the summary.
#' @param level Confidence level for individual-weight intervals and pointwise
#'   effect intervals. Curve intervals are not simultaneous confidence bands.
#' @param effects Effect estimand to include. `"rowwise"` (the default) computes
#'   the average observed-row predictive finite difference by changing the focal
#'   predictor within each observed covariate row and then averaging. `"partial"`
#'   computes the average partial-dependence finite-difference contrast used by
#'   [pce()], which additionally averages over the observed focal-predictor
#'   values. These estimands generally differ when the fitted response contains
#'   interactions. `"none"` omits effect summaries.
#' @name statnnet-methods
NULL

#' @rdname statnnet-methods
#' @export
print.statnnet <- function(x, ...) {
  cat("statnnet: statistical summary of a fitted nnet model\n")
  cat("Formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat(
    sprintf(
      "Architecture: %d-%d-%d; %d weights; n = %d\n",
      x$architecture$inputs,
      x$architecture$hidden,
      x$architecture$outputs,
      length(x$weights),
      nrow(x$x)
    )
  )
  cat(sprintf("Response: %s; decay: %g\n", x$response, x$decay))
  cat(sprintf("nnet convergence code: %d\n", x$diagnostics$convergence))
  if (isTRUE(x$diagnostics$covariance_available)) {
    cat(
      sprintf(
        "Covariance: available (reciprocal condition number %.3g)\n",
        x$diagnostics$information_rcond
      )
    )
  } else {
    cat("Covariance: unavailable - ", x$diagnostics$covariance_reason, "\n", sep = "")
  }
  invisible(x)
}

#' @rdname statnnet-methods
#' @export
coef.statnnet <- function(object, ...) {
  object$weights
}

#' @rdname statnnet-methods
#' @export
vcov.statnnet <- function(object, ...) {
  if (!isTRUE(object$diagnostics$covariance_available) || is.null(object$covariance)) {
    stop(
      paste0(
        "The covariance matrix is unavailable: ",
        object$diagnostics$covariance_reason,
        "."
      ),
      call. = FALSE
    )
  }
  object$covariance
}

#' @rdname statnnet-methods
#' @export
predict.statnnet <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(stats::predict(object$fit, ...))
  }
  stats::predict(object$fit, newdata = newdata, ...)
}

#' @rdname statnnet-methods
#' @export
anova.statnnet <- function(object, ...) {
  if (length(list(...)) > 0L) {
    stop("Model-comparison arguments are not supported by `anova.statnnet()`.", call. = FALSE)
  }
  rows <- lapply(object$groups, function(group) {
    indices <- group$weight_indices
    effective_df <- NA_real_
    statistic <- NA_real_
    p_value <- NA_real_
    status <- "covariance unavailable"

    if (isTRUE(object$diagnostics$covariance_available)) {
      group_covariance <- object$covariance[indices, indices, drop = FALSE]
      group_scale <- max(1, max(abs(group_covariance)))
      group_eigenvalues <- eigen(
        group_covariance,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      group_rcond <- base::rcond(group_covariance)
      effective_df <- sum(diag(object$shrinkage[indices, indices, drop = FALSE]))

      if (min(group_eigenvalues) <= object$covariance_tol * group_scale ||
          !is.finite(group_rcond) || group_rcond < object$covariance_tol) {
        status <- "group covariance is singular"
      } else if (!is.finite(effective_df) || effective_df <= 0) {
        status <- "effective degrees of freedom are not positive"
      } else {
        inverse_group <- chol2inv(chol(group_covariance))
        theta <- object$weights[indices]
        statistic <- drop(crossprod(theta, inverse_group %*% theta))
        p_value <- stats::pchisq(statistic, df = effective_df, lower.tail = FALSE)
        status <- "ok"
      }
    }

    data.frame(
      term = group$term,
      model_columns = paste(group$columns, collapse = ", "),
      n_weights = length(indices),
      effective_df = effective_df,
      wald_chisq = statistic,
      p_value = p_value,
      status = status,
      stringsAsFactors = FALSE
    )
  })
  result <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      term = character(),
      model_columns = character(),
      n_weights = integer(),
      effective_df = numeric(),
      wald_chisq = numeric(),
      p_value = numeric(),
      status = character(),
      stringsAsFactors = FALSE
    )
  }
  rownames(result) <- NULL
  attr(result, "heading") <- "Grouped input-to-hidden Wald summaries"
  class(result) <- c("anova.statnnet", "data.frame")
  result
}

#' @export
print.anova.statnnet <- function(x, ...) {
  cat(attr(x, "heading"), "\n", sep = "")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' @rdname statnnet-methods
#' @export
summary.statnnet <- function(object, weights = FALSE, level = 0.95,
                             effects = c("rowwise", "none", "partial"), ...) {
  if (length(list(...)) > 0L) {
    stop("Unused arguments were supplied to `summary.statnnet()`.", call. = FALSE)
  }
  if (!is.logical(weights) || length(weights) != 1L || is.na(weights)) {
    stop("`weights` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be between zero and one.", call. = FALSE)
  }
  effects <- match.arg(effects)

  uncertainty <- if (isTRUE(object$diagnostics$covariance_available)) "delta" else "none"
  pce_rows <- switch(
    effects,
    none = list(),
    rowwise = lapply(object$variables, function(variable) {
      .average_predictive_effect(
        object,
        variable = variable,
        uncertainty = uncertainty,
        level = level
      )
    }),
    partial = lapply(object$variables, function(variable) {
      pce(
        object,
        variable = variable,
        type = "average",
        uncertainty = uncertainty,
        level = level
      )
    })
  )
  pce_table <- if (length(pce_rows)) do.call(rbind, pce_rows) else data.frame()
  rownames(pce_table) <- NULL

  result <- list(
    call = object$call,
    formula = object$formula,
    response = object$response,
    architecture = object$architecture,
    objective = object$objective,
    diagnostics = object$diagnostics,
    effects = effects,
    pce = pce_table,
    grouped_wald = stats::anova(object),
    parameter_wald = if (weights) .parameter_wald(object, level = level) else NULL
  )
  class(result) <- "summary.statnnet"
  result
}

#' @export
print.summary.statnnet <- function(x, ...) {
  cat("statnnet summary\n")
  cat("Formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat(
    sprintf(
      "Response: %s; architecture: %d-%d-%d\n",
      x$response,
      x$architecture$inputs,
      x$architecture$hidden,
      x$architecture$outputs
    )
  )
  if (!isTRUE(x$diagnostics$covariance_available)) {
    cat("Uncertainty unavailable: ", x$diagnostics$covariance_reason, "\n", sep = "")
  }
  if (identical(x$effects, "rowwise")) {
    cat("\nAverage observed-row predictive finite-difference effects\n")
    print.data.frame(x$pce, row.names = FALSE)
  } else if (identical(x$effects, "partial")) {
    cat("\nAverage partial covariate effects\n")
    print.data.frame(x$pce, row.names = FALSE)
  } else {
    cat("\nEffect summaries: omitted\n")
  }
  cat("\n")
  print(x$grouped_wald)
  if (!is.null(x$parameter_wald)) {
    cat("\nIndividual-weight diagnostics\n")
    print.data.frame(x$parameter_wald, row.names = FALSE)
  }
  invisible(x)
}
