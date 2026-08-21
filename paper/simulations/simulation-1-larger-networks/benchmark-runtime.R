script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
simulation_dir <- dirname(script_path)

source(file.path(simulation_dir, "config.R"))
source(file.path(simulation_dir, "simulation-functions.R"))

arguments <- commandArgs(trailingOnly = TRUE)
value_of <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), arguments, value = TRUE)
  if (!length(hit)) return(default)
  as.integer(sub(paste0("^--", name, "="), "", hit[[1L]]))
}

benchmark_starts <- value_of("starts", 2L)
target_starts <- value_of("target-starts", 10L)
target_replicates <- value_of("target-replicates", 1000L)
config <- simulation_config(replicates = 1L, starts = benchmark_starts)
config$truth_n <- 100000L
config$evaluation_n <- 500L
validate_simulation_config(config)
scenarios <- simulation_scenarios(config)

if (!requireNamespace("nnet", quietly = TRUE) ||
    !requireNamespace("statnnet", quietly = TRUE)) {
  stop("Install nnet and the current statnnet source first.", call. = FALSE)
}

message("Preparing benchmark truth...")
truth <- prepare_simulation_truth(config)
rows <- vector("list", nrow(scenarios))

for (index in seq_len(nrow(scenarios))) {
  scenario <- scenarios[index, , drop = FALSE]
  message(sprintf(
    "Benchmarking %02d/%02d: %s, n=%d, q=%d",
    index, nrow(scenarios), scenario$outcome, scenario$n, scenario$q_fit
  ))
  result <- run_simulation_replication(config, scenario, 1L, truth)
  fit_per_start <- result$fit_elapsed_seconds / benchmark_starts
  estimated_per_replication <-
    fit_per_start * target_starts + result$analysis_elapsed_seconds
  rows[[index]] <- data.frame(
    scenario_id = scenario$scenario_id,
    outcome = scenario$outcome,
    n = scenario$n,
    q_fit = scenario$q_fit,
    benchmark_starts = benchmark_starts,
    benchmark_seconds = result$total_elapsed_seconds,
    fit_seconds_per_start = fit_per_start,
    analysis_seconds = result$analysis_elapsed_seconds,
    estimated_seconds_per_replication = estimated_per_replication,
    target_replicates = target_replicates,
    estimated_scenario_hours = estimated_per_replication * target_replicates / 3600,
    stringsAsFactors = FALSE
  )
}

benchmark <- do.call(rbind, rows)
output <- file.path(
  simulation_dir,
  paste0("runtime-benchmark-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
)
utils::write.csv(benchmark, output, row.names = FALSE)
total_hours <- sum(benchmark$estimated_scenario_hours)

print(benchmark, row.names = FALSE)
cat("\nEstimated full single-process runtime: ", format_duration(total_hours * 3600), "\n", sep = "")
cat("Idealised runtime with 4 scenario-level processes: ", format_duration(total_hours * 3600 / 4), "\n", sep = "")
cat("Idealised runtime with 8 scenario-level processes: ", format_duration(total_hours * 3600 / 8), "\n", sep = "")
cat("Benchmark details: ", output, "\n", sep = "")
