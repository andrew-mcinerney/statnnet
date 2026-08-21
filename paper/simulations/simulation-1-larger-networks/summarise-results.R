script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
simulation_dir <- dirname(script_path)

arguments <- commandArgs(trailingOnly = TRUE)
run_id_argument <- grep("^--run-id=", arguments, value = TRUE)
if (length(run_id_argument) != 1L) {
  stop("Supply exactly one --run-id=NAME argument.", call. = FALSE)
}
run_id <- sub("^--run-id=", "", run_id_argument)
run_dir <- file.path(simulation_dir, "results", run_id)
raw_files <- list.files(
  file.path(run_dir, "raw"),
  pattern = "^scenario-[0-9]+[.]rds$",
  full.names = TRUE
)
if (!length(raw_files)) stop("No scenario checkpoints were found.", call. = FALSE)

results <- do.call(rbind, lapply(raw_files, function(path) readRDS(path)$results))
results <- results[order(results$scenario_id, results$replicate), , drop = FALSE]

safe_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
safe_sd <- function(x) if (sum(is.finite(x)) > 1L) stats::sd(x[is.finite(x)]) else NA_real_
proportion_summary <- function(indicator) {
  indicator <- indicator[!is.na(indicator)]
  if (!length(indicator)) return(c(estimate = NA_real_, mc_se = NA_real_, n = 0))
  estimate <- mean(indicator)
  c(
    estimate = estimate,
    mc_se = sqrt(estimate * (1 - estimate) / length(indicator)),
    n = length(indicator)
  )
}

groups <- split(results, interaction(
  results$outcome,
  results$n,
  results$q_fit,
  drop = TRUE
))
summary_rows <- lapply(groups, function(data) {
  null_rejection <- ifelse(is.finite(data$null_wald_p), data$null_wald_p < 0.05, NA)
  active_rejection <- ifelse(is.finite(data$active_wald_p), data$active_wald_p < 0.05, NA)
  null_summary <- proportion_summary(null_rejection)
  power_summary <- proportion_summary(active_rejection)
  coverage_summary <- proportion_summary(data$pce_covered)
  empirical_pce_sd <- safe_sd(data$pce_estimate)
  mean_pce_se <- safe_mean(data$pce_std_error)

  data.frame(
    scenario_id = data$scenario_id[[1L]],
    outcome = data$outcome[[1L]],
    n = data$n[[1L]],
    q_fit = data$q_fit[[1L]],
    n_parameters = data$n_parameters[[1L]],
    n_per_parameter = data$n_per_parameter[[1L]],
    replicates = nrow(data),
    fit_success_rate = mean(data$fit_success),
    convergence_rate = mean(data$best_converged),
    covariance_availability_rate = mean(data$covariance_available),
    null_wald_reporting_rate = mean(is.finite(data$null_wald_p)),
    type1_error = null_summary[["estimate"]],
    type1_mc_se = null_summary[["mc_se"]],
    type1_n = null_summary[["n"]],
    active_wald_reporting_rate = mean(is.finite(data$active_wald_p)),
    power = power_summary[["estimate"]],
    power_mc_se = power_summary[["mc_se"]],
    power_n = power_summary[["n"]],
    pce_reporting_rate = mean(is.finite(data$pce_estimate)),
    pce_truth = safe_mean(data$pce_truth),
    pce_bias = safe_mean(data$pce_estimate - data$pce_truth),
    pce_empirical_sd = empirical_pce_sd,
    pce_mean_std_error = mean_pce_se,
    pce_se_ratio = mean_pce_se / empirical_pce_sd,
    pce_coverage = coverage_summary[["estimate"]],
    pce_coverage_mc_se = coverage_summary[["mc_se"]],
    pce_coverage_n = coverage_summary[["n"]],
    mean_prediction_rmse = safe_mean(data$prediction_rmse_to_mean),
    mean_elapsed_seconds = safe_mean(data$total_elapsed_seconds),
    unexpected_error_rate = mean(!is.na(data$error_message)),
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[order(
  summary_table$outcome,
  summary_table$n,
  summary_table$q_fit
), , drop = FALSE]
rownames(summary_table) <- NULL

utils::write.csv(results, file.path(run_dir, "combined-results.csv"), row.names = FALSE)
utils::write.csv(summary_table, file.path(run_dir, "simulation-summary.csv"), row.names = FALSE)

failure_labels <- ifelse(
  !is.na(results$error_message),
  paste(results$error_stage, results$error_message, sep = ": "),
  ifelse(
    !results$covariance_available & !is.na(results$covariance_reason),
    paste("covariance unavailable", results$covariance_reason, sep = ": "),
    "none"
  )
)
failure_table <- as.data.frame(table(failure_labels), stringsAsFactors = FALSE)
names(failure_table) <- c("result", "frequency")
failure_table <- failure_table[order(-failure_table$frequency), , drop = FALSE]
utils::write.csv(failure_table, file.path(run_dir, "failure-summary.csv"), row.names = FALSE)

print(summary_table, row.names = FALSE)
cat("\nWrote summaries to ", run_dir, "\n", sep = "")
