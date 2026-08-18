.validate_nnet_structure <- function(object) {
  if (isTRUE(object$softmax)) {
    stop("Softmax nnet models are not supported.", call. = FALSE)
  }
  if (isTRUE(object$censored)) {
    stop("Censored nnet models are not supported.", call. = FALSE)
  }
  if (is.null(object$n) || length(object$n) != 3L || any(object$n < 1L)) {
    stop("The fitted nnet architecture is invalid or cannot be reconstructed.", call. = FALSE)
  }
  if (object$n[3L] != 1L) {
    stop("Only nnet models with one output node are supported.", call. = FALSE)
  }

  expected_weights <- (object$n[1L] + 1L) * object$n[2L] +
    object$n[2L] + 1L
  if (length(object$wts) != expected_weights) {
    stop("Skip-layer connections and nonstandard nnet architectures are not supported.", call. = FALSE)
  }
  if (!is.null(object$call$mask)) {
    stop("nnet fits with an explicitly supplied `mask` are not supported.", call. = FALSE)
  }
  if (length(object$decay) != 1L || !is.finite(object$decay) || object$decay < 0) {
    stop("Only one non-negative scalar `decay` value is supported.", call. = FALSE)
  }
  if (any(!is.finite(object$wts))) {
    stop("The fitted nnet weights must all be finite.", call. = FALSE)
  }
  invisible(TRUE)
}

.nnet_response_type <- function(object) {
  if (isTRUE(object$entropy)) {
    if (object$nsunits != object$nunits) {
      stop("The entropy fit has an incompatible nnet output structure.", call. = FALSE)
    }
    return("binary")
  }
  if (object$nsunits == object$nunits - object$n[3L]) {
    return("continuous")
  }
  stop(
    "Only linear-output regression and single-output binary entropy fits are supported.",
    call. = FALSE
  )
}

.reconstruct_nnet_data <- function(object, formula, data, response,
                                   evaluation_environment) {
  data <- as.data.frame(data)
  supplied_terms <- stats::terms(formula, data = data)
  if (!is.null(object$terms)) {
    fitted_formula <- stats::formula(object$terms)
    supplied_formula <- stats::formula(supplied_terms)
    if (!identical(as.character(fitted_formula), as.character(supplied_formula))) {
      stop("`formula` does not match the formula stored in the fitted nnet object.", call. = FALSE)
    }
    terms <- object$terms
  } else {
    terms <- supplied_terms
  }
  model_frame_all <- stats::model.frame(
    terms,
    data = data,
    na.action = stats::na.pass,
    xlev = object$xlevels
  )
  complete <- stats::complete.cases(model_frame_all)
  model_frame <- model_frame_all[complete, , drop = FALSE]
  data_used <- data[complete, , drop = FALSE]

  n_fitted <- NROW(object$fitted.values)
  if (nrow(model_frame) != n_fitted) {
    stop(
      "The supplied formula and data do not reconstruct the observations used by nnet.",
      call. = FALSE
    )
  }
  fitted_rows <- rownames(object$fitted.values)
  if (!is.null(fitted_rows) && !identical(rownames(model_frame), fitted_rows)) {
    stop(
      "The rows of `data` are not in the same order as the observations used by nnet.",
      call. = FALSE
    )
  }

  if (!is.null(object$call$weights)) {
    case_weights <- tryCatch(
      eval(object$call$weights, envir = data, enclos = evaluation_environment),
      error = function(error) NULL
    )
    if (is.null(case_weights) || length(case_weights) != nrow(data)) {
      stop("The case weights used by nnet could not be reconstructed.", call. = FALSE)
    }
    case_weights <- case_weights[complete]
    if (any(!is.finite(case_weights)) || any(case_weights != 1)) {
      stop("Case weights other than one are not currently supported.", call. = FALSE)
    }
  }

  model_matrix_full <- stats::model.matrix(
    terms,
    data = model_frame,
    contrasts.arg = object$contrasts
  )
  intercept <- match("(Intercept)", colnames(model_matrix_full), nomatch = 0L)
  keep <- if (intercept > 0L) {
    seq_len(ncol(model_matrix_full)) != intercept
  } else {
    rep(TRUE, ncol(model_matrix_full))
  }
  x <- model_matrix_full[, keep, drop = FALSE]
  assign <- attr(model_matrix_full, "assign")[keep]

  if (ncol(x) != object$n[1L]) {
    stop(
      "The reconstructed model matrix has a different number of inputs from the fitted nnet object.",
      call. = FALSE
    )
  }
  if (!is.null(object$coefnames) && !identical(colnames(x), object$coefnames)) {
    stop(
      "The reconstructed model-matrix columns do not match the fitted nnet input order.",
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop("The reconstructed model matrix contains non-finite values.", call. = FALSE)
  }

  y_original <- stats::model.response(model_frame)
  response_levels <- NULL
  if (response == "binary") {
    if (is.factor(y_original)) {
      if (nlevels(y_original) != 2L) {
        stop("A binary response factor must have exactly two levels.", call. = FALSE)
      }
      response_levels <- levels(y_original)
      y <- as.numeric(y_original) - 1
    } else {
      y <- as.numeric(y_original)
      if (!all(y %in% c(0, 1))) {
        stop("A binary response must contain only 0 and 1.", call. = FALSE)
      }
      response_levels <- object$lev
    }
  } else {
    if (!is.numeric(y_original) || NCOL(y_original) != 1L) {
      stop("A continuous response must be one numeric column.", call. = FALSE)
    }
    y <- as.numeric(y_original)
  }
  if (any(!is.finite(y))) {
    stop("The response contains non-finite values.", call. = FALSE)
  }

  variables <- unique(all.vars(stats::delete.response(terms)))
  variables <- variables[variables %in% names(data_used)]

  list(
    data = data_used,
    terms = terms,
    contrasts = object$contrasts,
    xlevels = object$xlevels,
    x = x,
    y = y,
    response_levels = response_levels,
    assign = assign,
    term_labels = attr(terms, "term.labels"),
    column_names = colnames(x),
    variables = variables
  )
}

.sigmoid <- function(x) {
  stats::plogis(x)
}

.nn_components <- function(x, weights, hidden, response) {
  x <- as.matrix(x)
  n_inputs <- ncol(x)
  hidden_count <- (n_inputs + 1L) * hidden
  hidden_weights <- matrix(
    weights[seq_len(hidden_count)],
    nrow = hidden,
    ncol = n_inputs + 1L,
    byrow = TRUE
  )
  hidden_linear <- cbind(1, x) %*% t(hidden_weights)
  hidden_activation <- .sigmoid(hidden_linear)
  output_weights <- weights[hidden_count + seq_len(hidden + 1L)]
  output_linear <- drop(cbind(1, hidden_activation) %*% output_weights)
  fitted <- if (response == "binary") .sigmoid(output_linear) else output_linear

  list(
    fitted = fitted,
    hidden = hidden_activation,
    output_weights = output_weights
  )
}

.nn_predict_matrix <- function(x, weights, hidden, response) {
  .nn_components(x, weights, hidden, response)$fitted
}

.nn_gradient_matrix <- function(x, weights, hidden, response) {
  x <- as.matrix(x)
  components <- .nn_components(x, weights, hidden, response)
  n <- nrow(x)
  n_inputs <- ncol(x)
  n_parameters <- length(weights)
  derivative_output <- if (response == "binary") {
    components$fitted * (1 - components$fitted)
  } else {
    rep(1, n)
  }

  gradient <- matrix(0, nrow = n, ncol = n_parameters)
  x_with_bias <- cbind(1, x)
  for (node in seq_len(hidden)) {
    derivative_hidden <- derivative_output *
      components$output_weights[node + 1L] *
      components$hidden[, node] * (1 - components$hidden[, node])
    indices <- (node - 1L) * (n_inputs + 1L) + seq_len(n_inputs + 1L)
    gradient[, indices] <- derivative_hidden * x_with_bias
  }
  output_indices <- (n_inputs + 1L) * hidden + seq_len(hidden + 1L)
  gradient[, output_indices] <- derivative_output * cbind(1, components$hidden)
  gradient
}

.validate_nnet_objective <- function(object, y, fitted, response, tolerance) {
  penalty <- as.numeric(object$decay) * sum(as.numeric(object$wts)^2)
  if (response == "continuous") {
    unpenalised <- sum((y - fitted)^2)
    criterion <- unpenalised + penalty
    label <- "residual sum of squares"
  } else {
    probability <- pmin(pmax(fitted, .Machine$double.eps), 1 - .Machine$double.eps)
    unpenalised <- -sum(y * log(probability) + (1 - y) * log1p(-probability))
    criterion <- unpenalised + penalty
    label <- "negative Bernoulli log-likelihood"
  }
  difference <- criterion - as.numeric(object$value)
  agrees <- is.finite(difference) &&
    abs(difference) <= tolerance * max(1, abs(as.numeric(object$value)))
  if (!agrees) {
    stop(
      sprintf(
        "The reconstructed %s plus decay penalty does not agree with `object$value`.",
        label
      ),
      call. = FALSE
    )
  }
  list(
    stored = as.numeric(object$value),
    reconstructed = criterion,
    unpenalised = unpenalised,
    penalty = penalty,
    difference = difference,
    agrees = agrees,
    sigma2 = if (response == "continuous") unpenalised / length(y) else NULL
  )
}

.nnet_curvature <- function(object, x, y, response, rss, tolerance) {
  hessian_source <- if (is.null(object$Hessian)) "recomputed" else "stored"
  hessian_penalised <- if (is.null(object$Hessian)) {
    tryCatch(
      nnet::nnetHess(object, x, y),
      error = function(error) {
        stop(
          paste("The nnet Hessian could not be recomputed:", conditionMessage(error)),
          call. = FALSE
        )
      }
    )
  } else {
    object$Hessian
  }
  hessian_penalised <- as.matrix(hessian_penalised)
  n_parameters <- length(object$wts)
  if (!identical(dim(hessian_penalised), c(n_parameters, n_parameters))) {
    stop("The nnet Hessian has incompatible dimensions.", call. = FALSE)
  }
  if (any(!is.finite(hessian_penalised))) {
    stop("The nnet Hessian contains non-finite values.", call. = FALSE)
  }
  symmetry_error <- max(abs(hessian_penalised - t(hessian_penalised)))
  symmetry_scale <- max(1, max(abs(hessian_penalised)))
  if (symmetry_error > tolerance * symmetry_scale) {
    stop("The nnet Hessian is not symmetric within numerical tolerance.", call. = FALSE)
  }
  hessian_penalised <- (hessian_penalised + t(hessian_penalised)) / 2
  hessian_unpenalised <- hessian_penalised -
    2 * as.numeric(object$decay) * diag(n_parameters)

  if (response == "continuous") {
    sigma2 <- rss / nrow(x)
    if (!is.finite(sigma2) || sigma2 <= 0) {
      stop("The Gaussian residual variance estimate is not positive and finite.", call. = FALSE)
    }
    information_penalised <- hessian_penalised / (2 * sigma2)
    information_unpenalised <- hessian_unpenalised / (2 * sigma2)
  } else {
    information_penalised <- hessian_penalised
    information_unpenalised <- hessian_unpenalised
  }

  eigen_penalised <- eigen(information_penalised, symmetric = TRUE, only.values = TRUE)$values
  eigen_unpenalised <- eigen(information_unpenalised, symmetric = TRUE, only.values = TRUE)$values
  information_scale <- max(1, max(abs(information_penalised)))
  eigen_tolerance <- tolerance * information_scale
  information_rcond <- base::rcond(information_penalised)
  covariance <- NULL
  shrinkage <- NULL
  covariance_available <- TRUE
  covariance_reason <- NULL

  if (min(eigen_penalised) <= eigen_tolerance) {
    covariance_available <- FALSE
    covariance_reason <- "the penalised information matrix is not positive definite"
  } else if (!is.finite(information_rcond) || information_rcond < tolerance) {
    covariance_available <- FALSE
    covariance_reason <- "the penalised information matrix is numerically singular"
  } else if (min(eigen_unpenalised) < -eigen_tolerance) {
    covariance_available <- FALSE
    covariance_reason <- "the unpenalised information matrix has material negative curvature"
  }

  if (covariance_available) {
    inverse_penalised <- chol2inv(chol(information_penalised))
    covariance <- inverse_penalised %*% information_unpenalised %*% inverse_penalised
    covariance <- (covariance + t(covariance)) / 2
    shrinkage <- inverse_penalised %*% information_unpenalised
    covariance_scale <- max(1, max(abs(covariance)))
    if (any(diag(covariance) < -tolerance * covariance_scale)) {
      covariance_available <- FALSE
      covariance_reason <- "the covariance calculation produced negative variances"
      covariance <- NULL
      shrinkage <- NULL
    }
  }

  diagnostics <- list(
    hessian_source = hessian_source,
    hessian_finite = TRUE,
    hessian_symmetry_error = symmetry_error,
    hessian_penalised_min_eigenvalue = min(eigen_penalised),
    hessian_unpenalised_min_eigenvalue = min(eigen_unpenalised),
    information_rcond = information_rcond,
    conditioning_warning = is.finite(information_rcond) &&
      information_rcond < sqrt(tolerance),
    covariance_available = covariance_available,
    covariance_reason = covariance_reason
  )

  list(
    hessian_penalised = hessian_penalised,
    hessian_unpenalised = hessian_unpenalised,
    information_penalised = information_penalised,
    information_unpenalised = information_unpenalised,
    covariance = covariance,
    shrinkage = shrinkage,
    diagnostics = diagnostics
  )
}

.make_term_groups <- function(assign, term_labels, column_names,
                              n_inputs, n_hidden) {
  group_ids <- unique(assign)
  group_ids <- group_ids[group_ids > 0L]
  groups <- lapply(group_ids, function(group_id) {
    columns <- which(assign == group_id)
    weight_indices <- sort(unlist(lapply(columns, function(column) {
      (seq_len(n_hidden) - 1L) * (n_inputs + 1L) + column + 1L
    })))
    list(
      term = term_labels[group_id],
      columns = column_names[columns],
      column_indices = columns,
      weight_indices = weight_indices
    )
  })
  names(groups) <- vapply(groups, `[[`, character(1), "term")
  groups
}

.build_model_matrix <- function(object, newdata) {
  terms <- stats::delete.response(object$terms)
  model_frame <- stats::model.frame(
    terms,
    data = as.data.frame(newdata),
    na.action = stats::na.pass,
    xlev = object$xlevels
  )
  if (any(!stats::complete.cases(model_frame))) {
    stop("`newdata` contains missing values in model variables.", call. = FALSE)
  }
  model_matrix <- stats::model.matrix(
    terms,
    data = model_frame,
    contrasts.arg = object$contrasts
  )
  intercept <- match("(Intercept)", colnames(model_matrix), nomatch = 0L)
  if (intercept > 0L) {
    model_matrix <- model_matrix[, -intercept, drop = FALSE]
  }
  if (!identical(colnames(model_matrix), colnames(object$x))) {
    stop(
      "`newdata` does not reproduce the fitted model-matrix columns and ordering.",
      call. = FALSE
    )
  }
  if (any(!is.finite(model_matrix))) {
    stop("The new model matrix contains non-finite values.", call. = FALSE)
  }
  model_matrix
}

.parameter_wald <- function(object, level = 0.95) {
  estimate <- object$weights
  result <- data.frame(
    parameter = names(estimate),
    estimate = unname(estimate),
    std_error = NA_real_,
    wald_chisq = NA_real_,
    p_value = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(object$diagnostics$covariance_available)) {
    return(result)
  }
  variance <- diag(object$covariance)
  valid <- is.finite(variance) & variance > 0
  result$std_error[valid] <- sqrt(variance[valid])
  result$wald_chisq[valid] <- estimate[valid]^2 / variance[valid]
  result$p_value[valid] <- stats::pchisq(
    result$wald_chisq[valid],
    df = 1,
    lower.tail = FALSE
  )
  critical <- stats::qnorm(1 - (1 - level) / 2)
  result$conf_low[valid] <- estimate[valid] - critical * result$std_error[valid]
  result$conf_high[valid] <- estimate[valid] + critical * result$std_error[valid]
  result
}
