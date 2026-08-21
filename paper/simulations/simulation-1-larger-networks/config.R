simulation_config <- function(replicates = 1000L, starts = 10L) {
  list(
    simulation = "larger-and-misspecified-architectures",
    version = 1L,
    base_seed = 20260821L,
    replicates = as.integer(replicates),
    starts = as.integer(starts),
    p = 8L,
    active_predictors = paste0("x", 1:5),
    null_predictors = paste0("x", 6:8),
    active_wald_term = "x1",
    null_wald_term = "x8",
    pce_variable = "x1",
    pce_value = 0,
    pce_step = 1,
    q_true = 4L,
    q_fit = c(2L, 4L, 8L, 12L),
    sample_sizes = c(500L, 1000L, 2000L),
    outcomes = c("gaussian", "binary"),
    predictor_rho = 0.3,
    gaussian_noise_sd = 1,
    decay = 0.01,
    maxit = 1000L,
    abstol = 1e-8,
    reltol = 1e-8,
    initial_weight_range = 0.7,
    covariance_tol = 1e-10,
    objective_tol = 1e-6,
    truth_n = 1000000L,
    evaluation_n = 2000L,
    checkpoint_every = 10L,
    progress_every = 10L
  )
}

simulation_scenarios <- function(config) {
  grid <- expand.grid(
    outcome = config$outcomes,
    n = config$sample_sizes,
    q_fit = config$q_fit,
    stringsAsFactors = FALSE
  )
  grid$scenario_id <- seq_len(nrow(grid))
  grid$n_parameters <- grid$q_fit * (config$p + 2L) + 1L
  grid$n_per_parameter <- grid$n / grid$n_parameters
  grid[, c(
    "scenario_id", "outcome", "n", "q_fit",
    "n_parameters", "n_per_parameter"
  )]
}
