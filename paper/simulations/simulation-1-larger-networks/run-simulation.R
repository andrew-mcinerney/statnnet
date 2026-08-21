script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
simulation_dir <- dirname(script_path)
repository_dir <- normalizePath(file.path(simulation_dir, "..", "..", ".."))

source(file.path(simulation_dir, "config.R"))
source(file.path(simulation_dir, "simulation-functions.R"))

parse_arguments <- function(arguments) {
  defaults <- list(
    run_id = format(Sys.time(), "run-%Y%m%d-%H%M%S"),
    replicates = 1000L,
    starts = 10L,
    scenario_id = "all",
    resume = FALSE,
    overwrite = FALSE,
    checkpoint_every = 10L,
    progress_every = 10L,
    truth_n = 1000000L,
    evaluation_n = 2000L
  )
  for (argument in arguments) {
    if (argument == "--help") {
      cat(
        "Usage: Rscript run-simulation.R [options]\n\n",
        "--run-id=NAME             Results subdirectory name\n",
        "--replicates=N            Replicates per selected scenario (default 1000)\n",
        "--starts=N                Random starts per fitted network (default 10)\n",
        "--scenario-id=all|1,2     Scenario IDs to run\n",
        "--resume=true             Resume an existing run\n",
        "--overwrite=true          Back up and restart an existing run\n",
        "--checkpoint-every=N      Replicates between checkpoint writes\n",
        "--progress-every=N        Replicates between progress messages\n",
        "--truth-n=N               Reference sample for the true PCE\n",
        "--evaluation-n=N          Fixed test covariate sample size\n",
        sep = ""
      )
      quit(save = "no", status = 0L)
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L) stop("Invalid argument: ", argument, call. = FALSE)
    key <- gsub("-", "_", pieces[[1L]], fixed = TRUE)
    if (!key %in% names(defaults)) stop("Unknown argument: ", argument, call. = FALSE)
    value <- pieces[[2L]]
    if (is.logical(defaults[[key]])) {
      defaults[[key]] <- tolower(value) %in% c("true", "1", "yes")
    } else if (is.integer(defaults[[key]])) {
      defaults[[key]] <- as.integer(value)
    } else {
      defaults[[key]] <- value
    }
  }
  defaults
}

options <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("nnet", quietly = TRUE) ||
    !requireNamespace("statnnet", quietly = TRUE)) {
  stop(
    "Install the current statnnet source and nnet before running the simulation.",
    call. = FALSE
  )
}

config <- simulation_config(options$replicates, options$starts)
config$checkpoint_every <- options$checkpoint_every
config$progress_every <- options$progress_every
config$truth_n <- options$truth_n
config$evaluation_n <- options$evaluation_n
validate_simulation_config(config)
scenarios <- simulation_scenarios(config)
simulation_source_files <- file.path(
  simulation_dir,
  c("config.R", "simulation-functions.R", "run-simulation.R")
)
source_checksums <- data.frame(
  file = basename(simulation_source_files),
  md5 = unname(tools::md5sum(simulation_source_files)),
  stringsAsFactors = FALSE
)

selected_ids <- if (identical(options$scenario_id, "all")) {
  scenarios$scenario_id
} else {
  as.integer(strsplit(options$scenario_id, ",", fixed = TRUE)[[1L]])
}
if (anyNA(selected_ids) || !all(selected_ids %in% scenarios$scenario_id)) {
  stop("`--scenario-id` contains an invalid scenario number.", call. = FALSE)
}

results_root <- file.path(simulation_dir, "results")
run_dir <- file.path(results_root, options$run_id)
if (dir.exists(run_dir) && !options$resume && !options$overwrite) {
  stop(
    "The run directory already exists. Use --resume=true or a new --run-id.",
    call. = FALSE
  )
}
if (dir.exists(run_dir) && options$overwrite) {
  backup <- paste0(run_dir, "-backup-", format(Sys.time(), "%Y%m%d-%H%M%S"))
  if (!file.rename(run_dir, backup)) {
    stop("Could not back up the existing run directory.", call. = FALSE)
  }
  message("Existing results moved to ", backup)
}
dir.create(file.path(run_dir, "raw"), recursive = TRUE, showWarnings = FALSE)

config_file <- file.path(run_dir, "config.rds")
if (file.exists(config_file) && options$resume) {
  previous_config <- readRDS(config_file)
  fields <- c(
    "version", "base_seed", "starts", "p", "q_true", "q_fit",
    "sample_sizes", "outcomes", "predictor_rho", "decay", "maxit",
    "truth_n", "evaluation_n"
  )
  if (!identical(previous_config[fields], config[fields])) {
    stop("The requested settings do not match the saved run configuration.", call. = FALSE)
  }
  previous_checksums <- utils::read.csv(
    file.path(run_dir, "source-checksums.csv"),
    stringsAsFactors = FALSE
  )
  if (!identical(previous_checksums, source_checksums)) {
    stop(
      "The simulation source has changed since this run was created. Start a new run ID.",
      call. = FALSE
    )
  }
} else {
  atomic_save_rds(config, config_file)
  utils::write.csv(scenarios, file.path(run_dir, "scenario-grid.csv"), row.names = FALSE)
  utils::write.csv(
    source_checksums,
    file.path(run_dir, "source-checksums.csv"),
    row.names = FALSE
  )
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "session-info.txt"))
  git_commit <- tryCatch(
    system2(
      "git",
      c("-C", shQuote(repository_dir), "rev-parse", "HEAD"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(condition) "unavailable"
  )
  git_status <- tryCatch(
    system2(
      "git",
      c("-C", shQuote(repository_dir), "status", "--short"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(condition) "unavailable"
  )
  writeLines(
    c(
      paste("Repository:", repository_dir),
      paste("Git commit:", paste(git_commit, collapse = " ")),
      paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste("Command:", paste(commandArgs(), collapse = " "))
    ),
    file.path(run_dir, "run-metadata.txt")
  )
  writeLines(git_status, file.path(run_dir, "git-status.txt"))
}

truth_file <- file.path(run_dir, "truth.rds")
truth <- if (file.exists(truth_file)) {
  readRDS(truth_file)
} else {
  message("Computing fixed truth and prediction reference samples...")
  value <- prepare_simulation_truth(config)
  atomic_save_rds(value, truth_file)
  value
}
utils::write.csv(truth$pce, file.path(run_dir, "pce-truth.csv"), row.names = FALSE)

for (scenario_id in selected_ids) {
  scenario <- scenarios[scenarios$scenario_id == scenario_id, , drop = FALSE]
  checkpoint <- file.path(
    run_dir,
    "raw",
    sprintf("scenario-%02d.rds", scenario_id)
  )
  stored <- if (file.exists(checkpoint) && options$resume) {
    readRDS(checkpoint)$results
  } else {
    NULL
  }
  completed <- if (is.null(stored)) integer() else stored$replicate
  remaining <- setdiff(seq_len(config$replicates), completed)
  message(
    sprintf(
      "Scenario %02d/%02d: outcome=%s, n=%d, q_fit=%d; %d replicates remaining.",
      scenario_id,
      nrow(scenarios),
      scenario$outcome,
      scenario$n,
      scenario$q_fit,
      length(remaining)
    )
  )
  scenario_started <- proc.time()[["elapsed"]]

  for (index in seq_along(remaining)) {
    replicate_id <- remaining[[index]]
    result <- tryCatch(
      run_simulation_replication(config, scenario, replicate_id, truth),
      error = function(condition) {
        failed <- empty_replication_result(config, scenario, replicate_id)
        failed$error_stage <- "unexpected runner error"
        failed$error_message <- conditionMessage(condition)
        failed
      }
    )
    stored <- rbind(stored, result)

    should_checkpoint <- index %% config$checkpoint_every == 0L ||
      index == length(remaining)
    if (should_checkpoint) {
      stored <- stored[order(stored$replicate), , drop = FALSE]
      atomic_save_rds(list(scenario = scenario, results = stored), checkpoint)
    }

    should_report <- index %% config$progress_every == 0L ||
      index == length(remaining)
    if (should_report) {
      elapsed <- proc.time()[["elapsed"]] - scenario_started
      rate <- elapsed / index
      eta <- rate * (length(remaining) - index)
      message(
        sprintf(
          "  completed %d/%d remaining replicates; elapsed %s; scenario ETA %s",
          index,
          length(remaining),
          format_duration(elapsed),
          format_duration(eta)
        )
      )
    }
  }
}

marker <- file.path(
  run_dir,
  paste0("TASK-COMPLETE-scenarios-", paste(selected_ids, collapse = "-"), ".txt")
)
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), marker)
message("Requested scenarios completed. Results: ", run_dir)
