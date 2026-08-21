script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
simulation_dir <- dirname(script_path)

source(file.path(simulation_dir, "config.R"))
source(file.path(simulation_dir, "simulation-functions.R"))

if (!requireNamespace("nnet", quietly = TRUE) ||
    !requireNamespace("statnnet", quietly = TRUE)) {
  stop("Install nnet and the current statnnet source first.", call. = FALSE)
}

config <- simulation_config(replicates = 1L, starts = 1L)
config$truth_n <- 10000L
config$evaluation_n <- 200L
validate_simulation_config(config)
scenarios <- simulation_scenarios(config)
stopifnot(nrow(scenarios) == 24L)
stopifnot(all(scenarios$n_parameters == scenarios$q_fit * (config$p + 2L) + 1L))

parameters <- true_network_parameters()
stopifnot(all(parameters$hidden_weights[, 7:9, drop = FALSE] == 0))
stopifnot(any(parameters$hidden_weights[, 2:6, drop = FALSE] != 0))

truth <- prepare_simulation_truth(config)
stopifnot(
  all(is.finite(truth$pce$pce_truth)),
  all(is.finite(truth$pce$truth_mc_se)),
  all(truth$pce$truth_mc_se > 0)
)

gaussian <- scenarios[
  scenarios$outcome == "gaussian" & scenarios$n == 500L & scenarios$q_fit == 2L,
  ,
  drop = FALSE
]
binary <- scenarios[
  scenarios$outcome == "binary" & scenarios$n == 500L & scenarios$q_fit == 2L,
  ,
  drop = FALSE
]

gaussian_result <- run_simulation_replication(config, gaussian, 1L, truth)
binary_result <- run_simulation_replication(config, binary, 1L, truth)
stopifnot(
  gaussian_result$fit_success,
  binary_result$fit_success,
  gaussian_result$statnnet_success,
  binary_result$statnnet_success,
  is.finite(gaussian_result$prediction_rmse_to_mean),
  is.finite(binary_result$prediction_rmse_to_mean)
)

gaussian_repeat <- run_simulation_replication(config, gaussian, 1L, truth)
comparison_columns <- setdiff(
  names(gaussian_result),
  c("fit_elapsed_seconds", "analysis_elapsed_seconds", "total_elapsed_seconds")
)
stopifnot(identical(
  gaussian_result[, comparison_columns, drop = FALSE],
  gaussian_repeat[, comparison_columns, drop = FALSE]
))

cat("Simulation setup validation passed.\n")
cat("Gaussian PCE truth:", truth$pce$pce_truth[truth$pce$outcome == "gaussian"], "\n")
cat("Binary PCE truth:", truth$pce$pce_truth[truth$pce$outcome == "binary"], "\n")
