#' Partial covariate effects for a fitted statnnet model
#'
#' Computes finite differences of the partial-dependence function. These are
#' descriptive summaries of the fitted response surface, not causal effects.
#' Correlated predictors can cause evaluations at poorly supported covariate
#' combinations.
#'
#' @param object A `statnnet` object.
#' @param variable Name of an original predictor.
#' @param values Values at which to evaluate a continuous predictor. By
#'   default, an evenly spaced grid is chosen so that both each baseline value
#'   and its finite-difference comparison remain within the observed range.
#'   User-supplied values are used unchanged and may therefore extrapolate.
#' @param d Finite-difference step for a continuous predictor. The default is
#'   one observed standard deviation. Binary predictors are always compared
#'   from 0 to 1 and factors against their first level.
#' @param length_out Number of grid points for a continuous predictor.
#' @param type `"curve"` for evaluations over `values`, or `"average"` for a
#'   point summary averaged over observed predictor values. Average effects
#'   may compare values beyond the observed range near a boundary.
#' @param uncertainty One of `"delta"`, `"simulation"`, or `"none"`.
#' @param level Pointwise confidence level. Curve intervals are not
#'   simultaneous confidence bands.
#' @param nsim Number of simulation draws when `uncertainty = "simulation"`.
#' @param seed Optional simulation seed.
#' @param by Optional name of a second predictor at whose values the PCE is
#'   evaluated. Continuous predictors use their observed mean minus and plus
#'   one standard deviation. Binary predictors use both fitted levels.
#'
#' @return A data frame of class `"statnnet_pce"`.
#' @export
pce <- function(object, variable, values = NULL, d = NULL, length_out = 101L,
                type = c("curve", "average"),
                uncertainty = c("delta", "simulation", "none"),
                level = 0.95, nsim = 1000L, seed = NULL, by = NULL) {
  if (!inherits(object, "statnnet")) {
    stop("`object` must inherit from \"statnnet\".", call. = FALSE)
  }
  if (!is.character(variable) || length(variable) != 1L ||
      !variable %in% object$variables) {
    stop("`variable` must name one predictor used in the fitted formula.", call. = FALSE)
  }
  type <- match.arg(type)
  uncertainty <- match.arg(uncertainty)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be between zero and one.", call. = FALSE)
  }
  length_out <- as.integer(length_out)
  if (length(length_out) != 1L || is.na(length_out) || length_out < 2L) {
    stop("`length_out` must be an integer of at least two.", call. = FALSE)
  }

  scenarios <- .pce_scenarios(
    object = object,
    variable = variable,
    values = values,
    d = d,
    length_out = length_out,
    average = identical(type, "average")
  )
  by_spec <- .pce_by_spec(object, variable = variable, by = by)
  if (!is.null(by_spec)) {
    if (identical(type, "average")) {
      stop("`by` is currently supported only for `type = \"curve\"`.", call. = FALSE)
    }
    n_scenarios <- length(scenarios$items)
    scenarios$items <- unlist(lapply(seq_along(by_spec$values), function(index) {
      lapply(scenarios$items, function(item) {
        item$low <- .set_pce_condition(
          item$low, variable = by, value = by_spec$values[[index]], object = object
        )
        item$high <- .set_pce_condition(
          item$high, variable = by, value = by_spec$values[[index]], object = object
        )
        item
      })
    }), recursive = FALSE)
    scenarios$metadata <- scenarios$metadata[
      rep(seq_len(nrow(scenarios$metadata)), times = length(by_spec$values)),
      ,
      drop = FALSE
    ]
    scenarios$metadata$by_variable <- by
    scenarios$metadata$by_value <- rep(
      unlist(by_spec$values, use.names = FALSE),
      each = n_scenarios
    )
    scenarios$metadata$by_label <- rep(by_spec$labels, each = n_scenarios)
    rownames(scenarios$metadata) <- NULL
  }
  evaluated <- lapply(scenarios$items, function(item) {
    x_low <- .build_model_matrix(object, item$low)
    x_high <- .build_model_matrix(object, item$high)
    estimate <- mean(.nn_predict_matrix(
      x_high, object$weights, object$architecture$hidden, object$response
    )) - mean(.nn_predict_matrix(
      x_low, object$weights, object$architecture$hidden, object$response
    ))
    gradient <- colMeans(.nn_gradient_matrix(
      x_high, object$weights, object$architecture$hidden, object$response
    ) - .nn_gradient_matrix(
      x_low, object$weights, object$architecture$hidden, object$response
    ))
    list(estimate = estimate, gradient = gradient, x_low = x_low, x_high = x_high)
  })

  if (isTRUE(scenarios$collapse_average)) {
    estimate <- mean(vapply(evaluated, `[[`, numeric(1), "estimate"))
    gradient <- colMeans(do.call(rbind, lapply(evaluated, `[[`, "gradient")))
    evaluated <- list(list(estimate = estimate, gradient = gradient))
    metadata <- scenarios$average_metadata
  } else {
    metadata <- scenarios$metadata
  }

  result <- data.frame(
    variable = variable,
    contrast = metadata$contrast,
    value = metadata$value,
    step = metadata$step,
    estimate = vapply(evaluated, `[[`, numeric(1), "estimate"),
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(by_spec)) {
    result$by_variable <- metadata$by_variable
    result$by_value <- metadata$by_value
    result$by_label <- metadata$by_label
  }

  covariance_available <- isTRUE(object$diagnostics$covariance_available)
  if (uncertainty != "none" && !covariance_available) {
    warning(
      paste0(
        "PCE uncertainty is unavailable: ",
        object$diagnostics$covariance_reason,
        "."
      ),
      call. = FALSE
    )
  } else if (uncertainty == "delta") {
    gradients <- do.call(rbind, lapply(evaluated, `[[`, "gradient"))
    variances <- rowSums((gradients %*% object$covariance) * gradients)
    variance_scale <- max(1, max(abs(variances)))
    invalid <- variances < -object$covariance_tol * variance_scale
    if (any(invalid)) {
      stop("The delta-method calculation produced a negative variance.", call. = FALSE)
    }
    variances[variances < 0] <- 0
    result$std_error <- sqrt(variances)
    critical <- stats::qnorm(1 - (1 - level) / 2)
    result$conf_low <- result$estimate - critical * result$std_error
    result$conf_high <- result$estimate + critical * result$std_error
  } else if (uncertainty == "simulation") {
    nsim <- as.integer(nsim)
    if (length(nsim) != 1L || is.na(nsim) || nsim < 2L) {
      stop("`nsim` must be an integer of at least two.", call. = FALSE)
    }
    simulated_weights <- .simulate_weights(object, nsim = nsim, seed = seed)
    simulated_effects <- vapply(seq_len(nrow(result)), function(index) {
      if (isTRUE(scenarios$collapse_average)) {
        effects <- vapply(scenarios$items, function(item) {
          x_low <- .build_model_matrix(object, item$low)
          x_high <- .build_model_matrix(object, item$high)
          .pce_simulated_effects(object, simulated_weights, x_low, x_high)
        }, numeric(nsim))
        rowMeans(effects)
      } else {
        .pce_simulated_effects(
          object,
          simulated_weights,
          evaluated[[index]]$x_low,
          evaluated[[index]]$x_high
        )
      }
    }, numeric(nsim))
    if (is.null(dim(simulated_effects))) {
      simulated_effects <- matrix(simulated_effects, ncol = 1L)
    }
    alpha <- (1 - level) / 2
    result$std_error <- apply(simulated_effects, 2L, stats::sd)
    result$conf_low <- apply(simulated_effects, 2L, stats::quantile, probs = alpha)
    result$conf_high <- apply(
      simulated_effects,
      2L,
      stats::quantile,
      probs = 1 - alpha
    )
  }

  class(result) <- c("statnnet_pce", "data.frame")
  attr(result, "type") <- type
  attr(result, "uncertainty") <- if (covariance_available) uncertainty else "none"
  attr(result, "level") <- level
  attr(result, "by") <- by
  result
}

.pce_by_spec <- function(object, variable, by) {
  if (is.null(by)) {
    return(NULL)
  }
  if (!is.character(by) || length(by) != 1L || is.na(by) ||
      !by %in% object$variables) {
    stop("`by` must name one predictor used in the fitted formula.", call. = FALSE)
  }
  if (identical(by, variable)) {
    stop("`by` must be different from `variable`.", call. = FALSE)
  }

  predictor <- object$training_data[[by]]
  fitted_levels <- object$xlevels[[by]]
  is_categorical <- is.factor(predictor) || !is.null(fitted_levels) ||
    is.character(predictor)
  if (is_categorical) {
    levels_used <- if (!is.null(fitted_levels)) {
      fitted_levels
    } else {
      levels(factor(predictor))
    }
    if (length(levels_used) != 2L) {
      stop("A categorical `by` predictor must have exactly two fitted levels.", call. = FALSE)
    }
    return(list(
      values = as.list(levels_used),
      labels = paste0(by, " = ", levels_used)
    ))
  }

  if (is.logical(predictor)) {
    return(list(
      values = list(FALSE, TRUE),
      labels = paste0(by, " = ", c(0, 1))
    ))
  }
  if (!is.numeric(predictor) || any(!is.finite(predictor))) {
    stop("`by` must be a continuous or binary predictor.", call. = FALSE)
  }
  unique_values <- sort(unique(predictor))
  if (length(unique_values) == 2L && all(unique_values == c(0, 1))) {
    return(list(
      values = as.list(c(0, 1)),
      labels = paste0(by, " = ", c(0, 1))
    ))
  }

  predictor_sd <- stats::sd(predictor)
  if (!is.finite(predictor_sd) || predictor_sd == 0) {
    stop("A continuous `by` predictor must have a non-zero standard deviation.", call. = FALSE)
  }
  conditioning_values <- mean(predictor) + c(-1, 1) * predictor_sd
  list(
    values = as.list(conditioning_values),
    labels = sprintf(
      "%s = %.4g (%s1 SD)",
      by,
      conditioning_values,
      c("-", "+")
    )
  )
}

.set_pce_condition <- function(data, variable, value, object) {
  predictor <- object$training_data[[variable]]
  fitted_levels <- object$xlevels[[variable]]
  if (is.factor(predictor) || !is.null(fitted_levels) || is.character(predictor)) {
    levels_used <- if (!is.null(fitted_levels)) fitted_levels else levels(factor(predictor))
    data[[variable]] <- factor(
      rep(as.character(value), nrow(data)),
      levels = levels_used,
      ordered = is.ordered(predictor)
    )
  } else if (is.logical(predictor)) {
    data[[variable]] <- rep(as.logical(value), nrow(data))
  } else {
    data[[variable]] <- rep(as.numeric(value), nrow(data))
  }
  data
}

.average_predictive_effect <- function(object, variable, d = NULL,
                                       uncertainty = c("delta", "none"),
                                       level = 0.95) {
  uncertainty <- match.arg(uncertainty)
  data <- object$training_data
  predictor <- data[[variable]]
  factor_levels <- object$xlevels[[variable]]

  is_categorical <- is.factor(predictor) || !is.null(factor_levels) ||
    is.character(predictor)
  unique_values <- if (is.numeric(predictor)) sort(unique(predictor)) else NULL
  is_binary <- (is.logical(predictor)) ||
    (is.numeric(predictor) && all(is.finite(predictor)) &&
       length(unique_values) == 2L && all(unique_values == c(0, 1)))

  # Fixed categorical and binary contrasts already require only O(n) work and
  # have the same rowwise and partial-dependence interpretations.
  if (is_categorical || is_binary) {
    return(pce(
      object,
      variable = variable,
      type = "average",
      uncertainty = uncertainty,
      level = level
    ))
  }
  if (!is.numeric(predictor) || any(!is.finite(predictor))) {
    stop("Average effects require numeric, logical, or factor predictors.", call. = FALSE)
  }
  if (is.null(d)) {
    d <- stats::sd(predictor)
  }
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d == 0) {
    stop("`d` must be one non-zero finite number for a continuous predictor.", call. = FALSE)
  }

  low <- data
  high <- data
  high[[variable]] <- predictor + d
  x_low <- .build_model_matrix(object, low)
  x_high <- .build_model_matrix(object, high)
  estimate <- mean(
    .nn_predict_matrix(
      x_high, object$weights, object$architecture$hidden, object$response
    ) - .nn_predict_matrix(
      x_low, object$weights, object$architecture$hidden, object$response
    )
  )
  gradient <- colMeans(
    .nn_gradient_matrix(
      x_high, object$weights, object$architecture$hidden, object$response
    ) - .nn_gradient_matrix(
      x_low, object$weights, object$architecture$hidden, object$response
    )
  )

  result <- data.frame(
    variable = variable,
    contrast = "average observed-row finite difference",
    value = mean(predictor),
    step = d,
    estimate = estimate,
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    stringsAsFactors = FALSE
  )
  if (uncertainty == "delta") {
    variance <- drop(gradient %*% object$covariance %*% gradient)
    variance_scale <- max(1, abs(variance))
    if (variance < -object$covariance_tol * variance_scale) {
      stop("The delta-method calculation produced a negative variance.", call. = FALSE)
    }
    variance <- max(0, variance)
    result$std_error <- sqrt(variance)
    critical <- stats::qnorm(1 - (1 - level) / 2)
    result$conf_low <- result$estimate - critical * result$std_error
    result$conf_high <- result$estimate + critical * result$std_error
  }

  class(result) <- c("statnnet_pce", "data.frame")
  attr(result, "type") <- "average_rowwise"
  attr(result, "uncertainty") <- uncertainty
  attr(result, "level") <- level
  result
}

.pce_scenarios <- function(object, variable, values, d, length_out, average) {
  data <- object$training_data
  predictor <- data[[variable]]
  n <- nrow(data)
  factor_levels <- object$xlevels[[variable]]

  if (is.factor(predictor) || !is.null(factor_levels) || is.character(predictor)) {
    levels_used <- if (!is.null(factor_levels)) factor_levels else levels(factor(predictor))
    if (length(levels_used) < 2L) {
      stop("A categorical predictor must have at least two fitted levels.", call. = FALSE)
    }
    reference <- levels_used[1L]
    comparisons <- if (is.null(values)) levels_used[-1L] else as.character(values)
    if (any(!comparisons %in% levels_used) || any(comparisons == reference)) {
      stop("Categorical `values` must be non-reference fitted levels.", call. = FALSE)
    }
    items <- lapply(comparisons, function(comparison) {
      low <- data
      high <- data
      low[[variable]] <- factor(
        rep(reference, n),
        levels = levels_used,
        ordered = is.ordered(predictor)
      )
      high[[variable]] <- factor(
        rep(comparison, n),
        levels = levels_used,
        ordered = is.ordered(predictor)
      )
      list(low = low, high = high)
    })
    metadata <- data.frame(
      contrast = paste(comparisons, "-", reference),
      value = NA_real_,
      step = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(
      items = items,
      metadata = metadata,
      collapse_average = FALSE,
      average_metadata = NULL
    ))
  }

  if (is.logical(predictor)) {
    predictor <- as.numeric(predictor)
  }
  if (!is.numeric(predictor) || any(!is.finite(predictor))) {
    stop("PCEs currently support numeric, logical, and factor predictors.", call. = FALSE)
  }
  unique_values <- sort(unique(predictor))
  if (length(unique_values) == 2L && all(unique_values == c(0, 1))) {
    low <- data
    high <- data
    low[[variable]] <- 0
    high[[variable]] <- 1
    return(list(
      items = list(list(low = low, high = high)),
      metadata = data.frame(
        contrast = "1 - 0",
        value = 1,
        step = 1,
        stringsAsFactors = FALSE
      ),
      collapse_average = FALSE,
      average_metadata = NULL
    ))
  }

  if (is.null(d)) {
    d <- stats::sd(predictor)
  }
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d == 0) {
    stop("`d` must be one non-zero finite number for a continuous predictor.", call. = FALSE)
  }
  if (is.null(values)) {
    values <- if (average) {
      predictor
    } else {
      predictor_range <- range(predictor)
      grid_limits <- if (d > 0) {
        c(predictor_range[1L], predictor_range[2L] - d)
      } else {
        c(predictor_range[1L] - d, predictor_range[2L])
      }
      if (grid_limits[1L] >= grid_limits[2L]) {
        stop(
          paste(
            "The absolute finite-difference step must be smaller than the",
            "observed predictor range."
          ),
          call. = FALSE
        )
      }
      seq(grid_limits[1L], grid_limits[2L], length.out = length_out)
    }
  }
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop("Continuous `values` must be finite numbers.", call. = FALSE)
  }
  items <- lapply(values, function(value) {
    low <- data
    high <- data
    low[[variable]] <- value
    high[[variable]] <- value + d
    list(low = low, high = high)
  })
  metadata <- data.frame(
    contrast = NA_character_,
    value = values,
    step = d,
    stringsAsFactors = FALSE
  )
  list(
    items = items,
    metadata = metadata,
    collapse_average = average,
    average_metadata = data.frame(
      contrast = "average finite difference",
      value = mean(values),
      step = d,
      stringsAsFactors = FALSE
    )
  )
}

.simulate_weights <- function(object, nsim, seed) {
  covariance <- object$covariance
  eigen_covariance <- eigen(covariance, symmetric = TRUE)
  scale <- max(1, max(abs(covariance)))
  if (min(eigen_covariance$values) < -object$covariance_tol * scale) {
    stop("The covariance matrix is not positive semidefinite for simulation.", call. = FALSE)
  }
  eigenvalues <- eigen_covariance$values
  eigenvalues[eigenvalues < 0] <- 0
  if (!is.null(seed)) {
    set.seed(seed)
  }
  standard_normal <- matrix(
    stats::rnorm(nsim * length(object$weights)),
    nrow = nsim
  )
  covariance_factor <- diag(
    sqrt(eigenvalues),
    nrow = length(eigenvalues)
  ) %*% t(eigen_covariance$vectors)
  sweep(
    standard_normal %*% covariance_factor,
    2L,
    object$weights,
    "+"
  )
}

.pce_simulated_effects <- function(object, simulated_weights, x_low, x_high) {
  apply(simulated_weights, 1L, function(weights) {
    mean(.nn_predict_matrix(
      x_high, weights, object$architecture$hidden, object$response
    )) - mean(.nn_predict_matrix(
      x_low, weights, object$architecture$hidden, object$response
    ))
  })
}
