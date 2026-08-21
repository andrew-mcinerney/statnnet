validate_simulation_config <- function(config) {
  stopifnot(
    is.list(config),
    config$p == 8L,
    config$q_true == 4L,
    config$replicates >= 1L,
    config$starts >= 1L,
    config$truth_n >= 1000L,
    config$evaluation_n >= 100L,
    config$predictor_rho >= 0,
    config$predictor_rho < 1,
    config$decay >= 0
  )
  invisible(config)
}

seed_from_components <- function(base_seed, ...) {
  values <- as.numeric(unlist(list(...), use.names = FALSE))
  modulus <- 2147483646
  seed <- as.numeric(base_seed)
  for (value in values) {
    seed <- (seed * 104729 + value * 13007 + 7919) %% modulus
  }
  as.integer(seed + 1)
}

ar1_covariance <- function(p, rho) {
  outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))
}

generate_predictors <- function(n, config, seed, antithetic = FALSE) {
  set.seed(seed)
  if (antithetic) {
    half_n <- ceiling(n / 2)
    z <- matrix(stats::rnorm(half_n * config$p), nrow = half_n)
    z <- rbind(z, -z)[seq_len(n), , drop = FALSE]
  } else {
    z <- matrix(stats::rnorm(n * config$p), nrow = n)
  }
  x <- z %*% chol(ar1_covariance(config$p, config$predictor_rho))
  colnames(x) <- paste0("x", seq_len(config$p))
  x
}

true_network_parameters <- function() {
  hidden_weights <- rbind(
    c(-0.6, 0.8, -0.5, 0.3, 0.6, -0.4, 0, 0, 0),
    c(-0.2, 0.4, 0.7, -0.6, 0.2, 0.5, 0, 0, 0),
    c( 0.2, -0.3, 0.2, 0.8, -0.5, 0.6, 0, 0, 0),
    c( 0.6, 0.2, -0.4, 0.3, 0.7, -0.2, 0, 0, 0)
  )
  list(
    hidden_weights = hidden_weights,
    output_weights = c(0, 1.0, -0.8, 0.7, 0.9)
  )
}

raw_true_signal <- function(x) {
  parameters <- true_network_parameters()
  hidden_linear <- cbind(1, x) %*% t(parameters$hidden_weights)
  hidden_activation <- stats::plogis(hidden_linear)
  drop(cbind(1, hidden_activation) %*% parameters$output_weights)
}

calibrate_true_network <- function(config) {
  x <- generate_predictors(
    config$truth_n,
    config,
    seed = seed_from_components(config$base_seed, 11),
    antithetic = TRUE
  )
  signal <- raw_true_signal(x)
  list(
    raw_mean = mean(signal),
    raw_sd = stats::sd(signal)
  )
}

true_response_mean <- function(x, outcome, calibration) {
  latent <- (raw_true_signal(x) - calibration$raw_mean) / calibration$raw_sd
  if (identical(outcome, "gaussian")) latent else stats::plogis(latent)
}

prepare_simulation_truth <- function(config) {
  calibration <- calibrate_true_network(config)
  truth_x <- generate_predictors(
    config$truth_n,
    config,
    seed = seed_from_components(config$base_seed, 12),
    antithetic = TRUE
  )
  low <- high <- truth_x
  low[, config$pce_variable] <- config$pce_value
  high[, config$pce_variable] <- config$pce_value + config$pce_step

  pce <- lapply(config$outcomes, function(outcome) {
    differences <- true_response_mean(high, outcome, calibration) -
      true_response_mean(low, outcome, calibration)
    data.frame(
      outcome = outcome,
      pce_truth = mean(differences),
      truth_mc_se = stats::sd(differences) / sqrt(length(differences)),
      stringsAsFactors = FALSE
    )
  })

  evaluation_x <- generate_predictors(
    config$evaluation_n,
    config,
    seed = seed_from_components(config$base_seed, 13),
    antithetic = TRUE
  )

  list(
    calibration = calibration,
    pce = do.call(rbind, pce),
    evaluation_data = as.data.frame(evaluation_x),
    evaluation_mean = setNames(
      lapply(config$outcomes, function(outcome) {
        true_response_mean(evaluation_x, outcome, calibration)
      }),
      config$outcomes
    )
  )
}

training_data_seed <- function(config, outcome, n, replicate_id) {
  seed_from_components(
    config$base_seed,
    21,
    match(outcome, config$outcomes),
    match(as.integer(n), config$sample_sizes),
    as.integer(replicate_id)
  )
}

fit_start_seed <- function(config, scenario, replicate_id, start_id) {
  seed_from_components(
    config$base_seed,
    31,
    scenario$scenario_id,
    as.integer(replicate_id),
    as.integer(start_id)
  )
}

simulate_training_data <- function(config, truth, outcome, n, replicate_id) {
  set.seed(training_data_seed(config, outcome, n, replicate_id))
  z <- matrix(stats::rnorm(n * config$p), nrow = n)
  x <- z %*% chol(ar1_covariance(config$p, config$predictor_rho))
  colnames(x) <- paste0("x", seq_len(config$p))
  response_mean <- true_response_mean(x, outcome, truth$calibration)
  y <- if (identical(outcome, "gaussian")) {
    response_mean + stats::rnorm(n, sd = config$gaussian_noise_sd)
  } else {
    stats::rbinom(n, size = 1L, prob = response_mean)
  }
  data.frame(y = y, x, check.names = FALSE)
}

capture_conditions <- function(expression) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      expression,
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) condition
  )
  list(
    value = value,
    error = if (inherits(value, "error")) conditionMessage(value) else NA_character_,
    warnings = unique(warnings)
  )
}

fit_multistart_nnet <- function(data, config, scenario, replicate_id) {
  formula <- y ~ .
  fit_results <- vector("list", config$starts)
  fit_elapsed <- system.time({
    for (start_id in seq_len(config$starts)) {
      set.seed(fit_start_seed(config, scenario, replicate_id, start_id))
      arguments <- list(
        formula,
        data = data,
        size = scenario$q_fit,
        decay = config$decay,
        maxit = config$maxit,
        abstol = config$abstol,
        reltol = config$reltol,
        rang = config$initial_weight_range,
        Hess = FALSE,
        trace = FALSE,
        MaxNWts = max(1000L, scenario$n_parameters + 10L)
      )
      if (identical(scenario$outcome, "gaussian")) {
        arguments$linout <- TRUE
      } else {
        arguments$entropy <- TRUE
      }
      fit_results[[start_id]] <- capture_conditions(
        do.call(nnet::nnet, arguments)
      )
    }
  })[["elapsed"]]

  successful <- vapply(
    fit_results,
    function(result) inherits(result$value, "nnet") && is.finite(result$value$value),
    logical(1)
  )
  converged <- vapply(fit_results, function(result) {
    inherits(result$value, "nnet") &&
      is.finite(result$value$value) &&
      identical(as.integer(result$value$convergence), 0L)
  }, logical(1))

  candidates <- which(converged)
  if (!length(candidates)) candidates <- which(successful)
  if (!length(candidates)) {
    return(list(
      fit = NULL,
      fit_elapsed = fit_elapsed,
      best_start = NA_integer_,
      n_successful_starts = sum(successful),
      n_converged_starts = sum(converged),
      start_warning_count = sum(vapply(fit_results, function(x) length(x$warnings), integer(1))),
      fit_error = paste(unique(stats::na.omit(vapply(fit_results, `[[`, character(1), "error"))), collapse = " | ")
    ))
  }

  objectives <- vapply(
    fit_results[candidates],
    function(result) as.numeric(result$value$value),
    numeric(1)
  )
  best_start <- candidates[[which.min(objectives)]]
  list(
    fit = fit_results[[best_start]]$value,
    fit_elapsed = fit_elapsed,
    best_start = best_start,
    n_successful_starts = sum(successful),
    n_converged_starts = sum(converged),
    start_warning_count = sum(vapply(fit_results, function(x) length(x$warnings), integer(1))),
    fit_error = NA_character_
  )
}

wald_row <- function(table, term) {
  row <- table[table$term == term, , drop = FALSE]
  if (!nrow(row)) {
    return(list(p_value = NA_real_, effective_df = NA_real_, status = "term absent"))
  }
  list(
    p_value = row$p_value[[1L]],
    effective_df = row$effective_df[[1L]],
    status = row$status[[1L]]
  )
}

empty_replication_result <- function(config, scenario, replicate_id) {
  data.frame(
    scenario_id = scenario$scenario_id,
    replicate = as.integer(replicate_id),
    outcome = scenario$outcome,
    n = scenario$n,
    q_true = config$q_true,
    q_fit = scenario$q_fit,
    n_parameters = scenario$n_parameters,
    n_per_parameter = scenario$n_per_parameter,
    decay = config$decay,
    starts = config$starts,
    best_start = NA_integer_,
    n_successful_starts = 0L,
    n_converged_starts = 0L,
    start_warning_count = 0L,
    fit_success = FALSE,
    best_converged = FALSE,
    objective = NA_real_,
    weight_l2_norm = NA_real_,
    statnnet_success = FALSE,
    covariance_available = FALSE,
    covariance_reason = NA_character_,
    information_rcond = NA_real_,
    conditioning_warning = NA,
    active_wald_p = NA_real_,
    active_wald_df = NA_real_,
    active_wald_status = NA_character_,
    null_wald_p = NA_real_,
    null_wald_df = NA_real_,
    null_wald_status = NA_character_,
    pce_truth = NA_real_,
    pce_truth_mc_se = NA_real_,
    pce_estimate = NA_real_,
    pce_std_error = NA_real_,
    pce_conf_low = NA_real_,
    pce_conf_high = NA_real_,
    pce_covered = NA,
    prediction_rmse_to_mean = NA_real_,
    fit_elapsed_seconds = NA_real_,
    analysis_elapsed_seconds = NA_real_,
    total_elapsed_seconds = NA_real_,
    error_stage = NA_character_,
    error_message = NA_character_,
    statnnet_warnings = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_simulation_replication <- function(config, scenario, replicate_id, truth) {
  started <- proc.time()[["elapsed"]]
  result <- empty_replication_result(config, scenario, replicate_id)
  truth_row <- truth$pce[truth$pce$outcome == scenario$outcome, , drop = FALSE]
  result$pce_truth <- truth_row$pce_truth
  result$pce_truth_mc_se <- truth_row$truth_mc_se

  data <- simulate_training_data(
    config,
    truth,
    outcome = scenario$outcome,
    n = scenario$n,
    replicate_id = replicate_id
  )
  fitted <- fit_multistart_nnet(data, config, scenario, replicate_id)
  result$fit_elapsed_seconds <- fitted$fit_elapsed
  result$best_start <- fitted$best_start
  result$n_successful_starts <- fitted$n_successful_starts
  result$n_converged_starts <- fitted$n_converged_starts
  result$start_warning_count <- fitted$start_warning_count
  if (is.null(fitted$fit)) {
    result$error_stage <- "nnet fitting"
    result$error_message <- fitted$fit_error
    result$total_elapsed_seconds <- proc.time()[["elapsed"]] - started
    return(result)
  }

  result$fit_success <- TRUE
  result$best_converged <- identical(as.integer(fitted$fit$convergence), 0L)
  result$objective <- as.numeric(fitted$fit$value)
  result$weight_l2_norm <- sqrt(sum(as.numeric(fitted$fit$wts)^2))

  analysed <- capture_conditions(
    statnnet::statnnet(
      fitted$fit,
      formula = y ~ .,
      data = data,
      response = if (scenario$outcome == "gaussian") "continuous" else "binary",
      covariance_tol = config$covariance_tol,
      objective_tol = config$objective_tol
    )
  )
  if (inherits(analysed$value, "error")) {
    result$error_stage <- "statnnet augmentation"
    result$error_message <- analysed$error
    result$statnnet_warnings <- paste(analysed$warnings, collapse = " | ")
    result$total_elapsed_seconds <- proc.time()[["elapsed"]] - started
    result$analysis_elapsed_seconds <- result$total_elapsed_seconds - result$fit_elapsed_seconds
    return(result)
  }

  model <- analysed$value
  result$statnnet_success <- TRUE
  result$statnnet_warnings <- paste(analysed$warnings, collapse = " | ")
  result$covariance_available <- isTRUE(model$diagnostics$covariance_available)
  result$covariance_reason <- if (is.null(model$diagnostics$covariance_reason)) {
    NA_character_
  } else {
    model$diagnostics$covariance_reason
  }
  result$information_rcond <- model$diagnostics$information_rcond
  result$conditioning_warning <- isTRUE(model$diagnostics$conditioning_warning)

  wald <- stats::anova(model)
  active <- wald_row(wald, config$active_wald_term)
  null <- wald_row(wald, config$null_wald_term)
  result$active_wald_p <- active$p_value
  result$active_wald_df <- active$effective_df
  result$active_wald_status <- active$status
  result$null_wald_p <- null$p_value
  result$null_wald_df <- null$effective_df
  result$null_wald_status <- null$status

  if (result$covariance_available) {
    pce_result <- capture_conditions(
      statnnet::pce(
        model,
        variable = config$pce_variable,
        values = config$pce_value,
        d = config$pce_step,
        uncertainty = "delta"
      )
    )
    if (!inherits(pce_result$value, "error")) {
      result$pce_estimate <- pce_result$value$estimate[[1L]]
      result$pce_std_error <- pce_result$value$std_error[[1L]]
      result$pce_conf_low <- pce_result$value$conf_low[[1L]]
      result$pce_conf_high <- pce_result$value$conf_high[[1L]]
      result$pce_covered <- result$pce_conf_low <= result$pce_truth &&
        result$pce_conf_high >= result$pce_truth
    } else {
      result$error_stage <- "PCE calculation"
      result$error_message <- pce_result$error
    }
  }

  predicted <- as.numeric(stats::predict(
    model,
    newdata = truth$evaluation_data
  ))
  target <- truth$evaluation_mean[[scenario$outcome]]
  result$prediction_rmse_to_mean <- sqrt(mean((predicted - target)^2))
  result$total_elapsed_seconds <- proc.time()[["elapsed"]] - started
  result$analysis_elapsed_seconds <- result$total_elapsed_seconds - result$fit_elapsed_seconds
  result
}

atomic_save_rds <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    copied <- file.copy(temporary, path, overwrite = TRUE)
    unlink(temporary)
    if (!copied) stop("Could not update checkpoint: ", path, call. = FALSE)
  }
  invisible(path)
}

format_duration <- function(seconds) {
  if (!is.finite(seconds)) return("unknown")
  if (seconds < 60) return(sprintf("%.1f seconds", seconds))
  if (seconds < 3600) return(sprintf("%.1f minutes", seconds / 60))
  if (seconds < 86400) return(sprintf("%.1f hours", seconds / 3600))
  sprintf("%.1f days", seconds / 86400)
}
