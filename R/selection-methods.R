#' Methods for statnnet model-selection results
#'
#' @param x,object A `statnnet_selection` object.
#' @param newdata Optional new predictor data. For a formula selection this is
#'   a data frame; for a matrix selection it may be a matrix or data frame.
#' @param ... Additional arguments passed to the underlying method where
#'   documented.
#' @examples
#' set.seed(1)
#' example_data <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' example_data$y <- 1 + 2 * example_data$x1 + rnorm(60, sd = 0.4)
#' selection <- select_nnet(
#'   y ~ x1 + x2,
#'   data = example_data,
#'   max_size = 1,
#'   n_starts = 2,
#'   maxit = 500,
#'   seed = 2,
#'   progress = FALSE
#' )
#'
#' selection
#' summary(selection)
#' coef(selection)
#' head(predict(selection))
#' plot(selection)
#' @name statnnet-selection-methods
NULL

#' @rdname statnnet-selection-methods
#' @export
print.statnnet_selection <- function(x, ...) {
  cat("statnnet model selection\n")
  cat("Call:\n")
  print(x$call)
  cat(
    sprintf(
      "Selected model: %d input column%s, %d hidden unit%s, %d network weights\n",
      length(x$selected_columns),
      if (length(x$selected_columns) == 1L) "" else "s",
      x$size,
      if (x$size == 1L) "" else "s",
      length(x$fit$wts)
    )
  )
  cat(
    sprintf(
      "%s: %.3f; %d candidate models; %d converged starts\n",
      x$criterion,
      x$criterion_value,
      x$diagnostics$n_candidates,
      x$diagnostics$n_converged_starts
    )
  )
  cat("Selected terms: ", paste(x$selected_terms, collapse = ", "), "\n", sep = "")
  if (length(x$excluded_terms)) {
    cat("Excluded terms: ", paste(x$excluded_terms, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname statnnet-selection-methods
#' @export
summary.statnnet_selection <- function(object, ...) {
  if (length(list(...)) > 0L) {
    stop(
      "Unused arguments were supplied to `summary.statnnet_selection()`.",
      call. = FALSE
    )
  }
  term_table <- data.frame(
    term = c(object$selected_terms, object$excluded_terms),
    selected = c(
      rep(TRUE, length(object$selected_terms)),
      rep(FALSE, length(object$excluded_terms))
    ),
    stringsAsFactors = FALSE
  )
  result <- list(
    call = object$call,
    formula = object$formula,
    criterion = object$criterion,
    criterion_value = object$criterion_value,
    size = object$size,
    n_inputs = length(object$selected_columns),
    n_weights = length(object$fit$wts),
    terms = term_table,
    search_path = object$search_path,
    diagnostics = object$diagnostics
  )
  class(result) <- "summary.statnnet_selection"
  result
}

#' @export
print.summary.statnnet_selection <- function(x, ...) {
  cat("statnnet model-selection summary\n")
  cat("Selected formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat(
    sprintf(
      "Architecture: %d-%d-1; %d weights\n",
      x$n_inputs,
      x$size,
      x$n_weights
    )
  )
  cat(sprintf("%s: %.3f\n", x$criterion, x$criterion_value))
  cat(
    sprintf(
      "Search: %d candidates; %d failed candidates; %d failed starts\n",
      x$diagnostics$n_candidates,
      x$diagnostics$n_failed_candidates,
      x$diagnostics$n_failed_starts
    )
  )
  cat("\nTerms\n")
  print.data.frame(x$terms, row.names = FALSE)
  invisible(x)
}

#' @rdname statnnet-selection-methods
#' @export
coef.statnnet_selection <- function(object, ...) {
  stats::coef(object$fit, ...)
}

#' @rdname statnnet-selection-methods
#' @export
predict.statnnet_selection <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(stats::predict(object$fit, ...))
  }
  x <- .selection_newdata_matrix(object, newdata)
  stats::predict(object$fit, newdata = x, ...)
}

#' @rdname statnnet-selection-methods
#' @export
plot.statnnet_selection <- function(x, ...) {
  path <- selection_path(x)
  phase <- factor(path$phase, levels = unique(path$phase))
  colours <- seq_along(levels(phase))
  finite <- is.finite(path$criterion_value)
  if (!any(finite)) {
    stop("The selection path contains no finite criterion values.", call. = FALSE)
  }
  graphics::plot(
    path$evaluation[finite],
    path$criterion_value[finite],
    type = "b",
    pch = 19,
    col = colours[phase[finite]],
    xlab = "Candidate evaluation",
    ylab = x$criterion,
    ...
  )
  selected <- which(path$selected & finite)
  if (length(selected)) {
    graphics::points(
      path$evaluation[selected],
      path$criterion_value[selected],
      pch = 23,
      bg = "gold",
      cex = 1.5
    )
  }
  graphics::legend(
    "topright",
    legend = levels(phase),
    col = colours,
    pch = 19,
    bty = "n"
  )
  invisible(path)
}

#' Extract the model-selection path
#'
#' @param object A `statnnet_selection` object.
#' @return A data frame with one row per candidate model evaluated during the
#'   selection procedure.
#' @examples
#' set.seed(1)
#' example_data <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' example_data$y <- 1 + 2 * example_data$x1 + rnorm(60, sd = 0.4)
#' selection <- select_nnet(
#'   y ~ x1 + x2,
#'   data = example_data,
#'   max_size = 1,
#'   n_starts = 2,
#'   maxit = 500,
#'   seed = 2,
#'   progress = FALSE
#' )
#' selection_path(selection)
#' @export
selection_path <- function(object) {
  if (!inherits(object, "statnnet_selection")) {
    stop("`object` must inherit from \"statnnet_selection\".", call. = FALSE)
  }
  object$search_path
}

.selection_newdata_matrix <- function(object, newdata) {
  if (identical(object$interface, "formula")) {
    frame <- stats::model.frame(
      stats::delete.response(object$terms),
      data = as.data.frame(newdata),
      na.action = stats::na.pass,
      xlev = object$xlevels
    )
    matrix_full <- stats::model.matrix(
      stats::delete.response(object$terms),
      data = frame,
      contrasts.arg = object$contrasts
    )
    intercept <- match("(Intercept)", colnames(matrix_full), nomatch = 0L)
    keep <- if (intercept) seq_len(ncol(matrix_full)) != intercept else rep(TRUE, ncol(matrix_full))
    result <- matrix_full[, keep, drop = FALSE]
  } else {
    result <- as.matrix(newdata)
    if (!is.numeric(result)) {
      stop("Matrix-interface `newdata` must be numeric.", call. = FALSE)
    }
    if (!is.null(colnames(result)) &&
        all(object$selected_columns %in% colnames(result))) {
      result <- result[, object$selected_columns, drop = FALSE]
    } else if (ncol(result) == length(object$full_columns)) {
      match_index <- match(object$selected_columns, object$full_columns)
      result <- result[, match_index, drop = FALSE]
    } else if (ncol(result) != length(object$selected_columns)) {
      stop(
        "`newdata` does not contain the selected input columns.",
        call. = FALSE
      )
    }
  }
  if (ncol(result) != object$fit$n[1L]) {
    stop("`newdata` does not reconstruct the selected model inputs.", call. = FALSE)
  }
  result
}
