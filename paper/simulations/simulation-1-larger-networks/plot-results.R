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
summary_file <- file.path(run_dir, "simulation-summary.csv")
if (!file.exists(summary_file)) {
  stop("Run summarise-results.R before plotting.", call. = FALSE)
}

simulation_summary <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
required <- expand.grid(
  outcome = c("gaussian", "binary"),
  n = c(500L, 1000L, 2000L),
  q_fit = c(2L, 4L, 8L, 12L),
  stringsAsFactors = FALSE
)
if (nrow(merge(required, simulation_summary)) != nrow(required)) {
  warning("The diagnostic figure is based on an incomplete scenario grid.", call. = FALSE)
}

output_file <- file.path(run_dir, "simulation-diagnostics.pdf")
grDevices::pdf(output_file, width = 10, height = 6.5)
graphics::par(mfrow = c(2, 3), mar = c(3.7, 3.8, 2.4, 0.8), oma = c(0, 0, 1, 0))

metrics <- list(
  list(column = "type1_error", title = "Null Wald rejection", reference = 0.05),
  list(column = "pce_coverage", title = "PCE interval coverage", reference = 0.95),
  list(column = "covariance_availability_rate", title = "Covariance availability", reference = NA_real_)
)
colours <- c("#0072B2", "#D55E00", "#009E73")
point_types <- c(16, 17, 15)

for (outcome in c("gaussian", "binary")) {
  outcome_data <- simulation_summary[simulation_summary$outcome == outcome, , drop = FALSE]
  for (metric in metrics) {
    values <- outcome_data[[metric$column]]
    y_limits <- if (metric$column == "type1_error") {
      c(0, max(0.15, values, na.rm = TRUE))
    } else {
      c(0, 1)
    }
    graphics::plot(
      NA,
      xlim = range(required$q_fit),
      ylim = y_limits,
      xlab = "Fitted hidden nodes",
      ylab = metric$title,
      main = paste(tools::toTitleCase(outcome), metric$title),
      xaxt = "n"
    )
    graphics::axis(1, at = sort(unique(required$q_fit)))
    if (is.finite(metric$reference)) {
      graphics::abline(h = metric$reference, lty = 2, col = "grey40")
    }
    for (index in seq_along(colours)) {
      sample_size <- sort(unique(required$n))[[index]]
      line_data <- outcome_data[outcome_data$n == sample_size, , drop = FALSE]
      line_data <- line_data[order(line_data$q_fit), , drop = FALSE]
      graphics::lines(
        line_data$q_fit,
        line_data[[metric$column]],
        type = "b",
        col = colours[[index]],
        pch = point_types[[index]],
        lwd = 1.5
      )
    }
    if (identical(metric$column, "type1_error")) {
      graphics::legend(
        "topleft",
        legend = paste("n =", sort(unique(required$n))),
        col = colours,
        pch = point_types,
        lty = 1,
        bty = "n",
        cex = 0.8
      )
    }
  }
}

grDevices::dev.off()
cat("Wrote ", output_file, "\n", sep = "")
