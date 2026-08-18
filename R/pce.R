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
#'   default, an evenly spaced grid over the observed range is used.
#' @param d Finite-difference step for a continuous predictor. The default is
#'   one observed standard deviation. Binary predictors are always compared
#'   from 0 to 1 and factors against their first level.
#' @param length_out Number of grid points for a continuous predictor.
#' @param type `"curve"` for evaluations over `values`, or `"average"` for a
#'   point summary averaged over observed predictor values.
#' @param uncertainty One of `"delta"`, `"simulation"`, or `"none"`.
#' @param level Confidence level.
#' @param nsim Number of simulation draws when `uncertainty = "simulation"`.
#' @param seed Optional simulation seed.
#'
#' @return A data frame of class `"statnnet_pce"`.
#' @export
pce <- function(object, variable, values = NULL, d = NULL, length_out = 101L,
                type = c("curve", "average"),
                uncertainty = c("delta", "simulation", "none"),
                level = 0.95, nsim = 1000L, seed = NULL) {
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
      seq(min(predictor), max(predictor), length.out = length_out)
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
